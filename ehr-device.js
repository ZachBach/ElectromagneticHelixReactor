// EHRDevice — classifies the viewport/input as mobile · tablet · desktop and
// stays current across resize/orientation change/OS on-screen-keyboard events.
// Combines viewport size with pointer/touch capability rather than raw
// user-agent sniffing, so foldables, touch laptops, and desktop dev-tools
// device emulation all classify sensibly.
(function () {
  function classify() {
    const w = window.innerWidth, h = window.innerHeight;
    const touch = ('ontouchstart' in window) || navigator.maxTouchPoints > 0;
    const coarse = typeof matchMedia === 'function' && matchMedia('(pointer: coarse)').matches;
    const touchLike = touch || coarse;

    // Width drives the classification — laptop screens are routinely
    // <1024px *tall* (e.g. 1440x900), so height must never demote a
    // desktop to "tablet". Touch only breaks ties in the iPad-landscape
    // zone (1024-1280) where a non-touch desktop window of the same
    // width should stay "desktop".
    let type;
    if (w < 640) type = 'mobile';
    else if (w < 1024) type = 'tablet';
    else if (touchLike && w < 1280) type = 'tablet';
    else type = 'desktop';

    const ua = navigator.userAgent || '';
    let os = 'desktop';
    if (/iPad|iPhone|iPod/.test(ua)) os = 'iOS';
    else if (/Android/.test(ua)) os = 'Android';
    else if (/Win/.test(ua)) os = 'Windows';
    else if (/Mac/.test(ua)) os = 'macOS';
    else if (/Linux/.test(ua)) os = 'Linux';

    return {
      type,                              // 'mobile' | 'tablet' | 'desktop'
      touch: !!(touch || coarse),
      os,
      width: w,
      height: h,
      orientation: w >= h ? 'landscape' : 'portrait',
    };
  }

  let current = classify();
  const subs = new Set();
  let raf = null;

  function recompute() {
    raf = null;
    const next = classify();
    const changed = next.type !== current.type || next.width !== current.width ||
      next.height !== current.height || next.orientation !== current.orientation;
    current = next;
    if (changed) subs.forEach((fn) => fn(current));
  }
  function scheduleRecompute() {
    if (raf) return;
    raf = requestAnimationFrame(recompute);
  }

  window.addEventListener('resize', scheduleRecompute);
  window.addEventListener('orientationchange', scheduleRecompute);

  window.EHRDevice = {
    get() { return current; },
    // subscribe(fn) — fn is called with the new device object whenever
    // type/size/orientation actually changes. Returns an unsubscribe fn.
    subscribe(fn) {
      subs.add(fn);
      return () => subs.delete(fn);
    },
  };
})();
