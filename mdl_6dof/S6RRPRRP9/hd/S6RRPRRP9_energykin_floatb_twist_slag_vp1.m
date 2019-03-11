% Calculate kinetic energy for
% S6RRPRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRP9_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP9_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP9_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRRP9_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP9_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP9_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP9_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP9_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:27:24
% EndTime: 2019-03-09 12:27:29
% DurationCPUTime: 4.28s
% Computational Cost: add. (3457->391), mult. (5858->549), div. (0->0), fcn. (6961->12), ass. (0->179)
t380 = Icges(6,1) + Icges(7,1);
t379 = Icges(6,4) + Icges(7,4);
t378 = Icges(6,5) + Icges(7,5);
t377 = Icges(6,2) + Icges(7,2);
t376 = Icges(6,6) + Icges(7,6);
t375 = Icges(6,3) + Icges(7,3);
t374 = rSges(7,3) + qJ(6);
t305 = cos(pkin(6));
t310 = sin(qJ(1));
t312 = cos(qJ(2));
t341 = t310 * t312;
t309 = sin(qJ(2));
t313 = cos(qJ(1));
t342 = t309 * t313;
t269 = t305 * t342 + t341;
t334 = pkin(11) + qJ(4);
t298 = sin(t334);
t327 = cos(t334);
t303 = sin(pkin(6));
t344 = t303 * t313;
t241 = t269 * t327 - t298 * t344;
t340 = t312 * t313;
t343 = t309 * t310;
t268 = -t305 * t340 + t343;
t308 = sin(qJ(5));
t311 = cos(qJ(5));
t211 = -t241 * t308 + t268 * t311;
t350 = t268 * t308;
t212 = t241 * t311 + t350;
t326 = t303 * t327;
t240 = t269 * t298 + t313 * t326;
t373 = t376 * t211 + t378 * t212 + t375 * t240;
t271 = -t305 * t343 + t340;
t346 = t303 * t310;
t243 = t271 * t327 + t298 * t346;
t270 = t305 * t341 + t342;
t213 = -t243 * t308 + t270 * t311;
t349 = t270 * t308;
t214 = t243 * t311 + t349;
t242 = t271 * t298 - t310 * t326;
t372 = t376 * t213 + t378 * t214 + t375 * t242;
t371 = t377 * t211 + t379 * t212 + t376 * t240;
t370 = t377 * t213 + t379 * t214 + t376 * t242;
t369 = t379 * t211 + t380 * t212 + t378 * t240;
t368 = t379 * t213 + t380 * t214 + t378 * t242;
t258 = t305 * t298 + t309 * t326;
t345 = t303 * t312;
t238 = -t258 * t308 - t311 * t345;
t328 = t308 * t345;
t239 = t258 * t311 - t328;
t347 = t303 * t309;
t257 = t298 * t347 - t305 * t327;
t367 = t376 * t238 + t378 * t239 + t375 * t257;
t366 = t377 * t238 + t379 * t239 + t376 * t257;
t365 = t379 * t238 + t380 * t239 + t378 * t257;
t302 = sin(pkin(11));
t304 = cos(pkin(11));
t244 = -t269 * t302 - t304 * t344;
t329 = t302 * t344;
t245 = t269 * t304 - t329;
t197 = Icges(4,5) * t245 + Icges(4,6) * t244 + Icges(4,3) * t268;
t229 = Icges(3,4) * t269 - Icges(3,2) * t268 - Icges(3,6) * t344;
t364 = t197 - t229;
t246 = -t271 * t302 + t304 * t346;
t330 = t302 * t346;
t247 = t271 * t304 + t330;
t198 = Icges(4,5) * t247 + Icges(4,6) * t246 + Icges(4,3) * t270;
t230 = Icges(3,4) * t271 - Icges(3,2) * t270 + Icges(3,6) * t346;
t363 = t198 - t230;
t266 = -t302 * t347 + t304 * t305;
t348 = t302 * t305;
t267 = t304 * t347 + t348;
t221 = Icges(4,5) * t267 + Icges(4,6) * t266 - Icges(4,3) * t345;
t255 = Icges(3,6) * t305 + (Icges(3,4) * t309 + Icges(3,2) * t312) * t303;
t362 = t221 - t255;
t355 = pkin(8) * t305;
t354 = pkin(3) * t304;
t353 = pkin(5) * t311;
t351 = Icges(2,4) * t310;
t338 = rSges(7,1) * t212 + rSges(7,2) * t211 + pkin(5) * t350 + t374 * t240 + t241 * t353;
t337 = rSges(7,1) * t214 + rSges(7,2) * t213 + pkin(5) * t349 + t374 * t242 + t243 * t353;
t336 = rSges(7,1) * t239 + rSges(7,2) * t238 - pkin(5) * t328 + t374 * t257 + t258 * t353;
t335 = qJD(2) * t303;
t333 = V_base(5) * pkin(7) + V_base(1);
t281 = t310 * t335 + V_base(4);
t299 = V_base(6) + qJD(1);
t249 = qJD(4) * t270 + t281;
t282 = qJD(2) * t305 + t299;
t280 = -t313 * t335 + V_base(5);
t274 = t310 * pkin(1) - pkin(8) * t344;
t325 = -t274 * t299 + V_base(5) * t355 + t333;
t275 = pkin(1) * t313 + pkin(8) * t346;
t324 = V_base(4) * t274 - t275 * V_base(5) + V_base(3);
t248 = qJD(4) * t268 + t280;
t264 = -qJD(4) * t345 + t282;
t272 = (pkin(2) * t309 - qJ(3) * t312) * t303;
t323 = qJD(3) * t270 + t280 * t272 + t325;
t322 = t299 * t275 + V_base(2) + (-pkin(7) - t355) * V_base(4);
t237 = pkin(2) * t271 + qJ(3) * t270;
t321 = qJD(3) * t268 + t282 * t237 + t322;
t236 = t269 * pkin(2) + t268 * qJ(3);
t320 = -qJD(3) * t345 + t281 * t236 + t324;
t195 = -pkin(3) * t329 + pkin(9) * t268 + t269 * t354;
t226 = pkin(3) * t348 + (-pkin(9) * t312 + t309 * t354) * t303;
t319 = t280 * t226 + (-t195 - t236) * t282 + t323;
t196 = pkin(3) * t330 + pkin(9) * t270 + t271 * t354;
t318 = t282 * t196 + (-t226 - t272) * t281 + t321;
t317 = t281 * t195 + (-t196 - t237) * t280 + t320;
t206 = pkin(4) * t241 + pkin(10) * t240;
t225 = pkin(4) * t258 + pkin(10) * t257;
t316 = -t206 * t264 + t248 * t225 + t319;
t207 = pkin(4) * t243 + pkin(10) * t242;
t315 = t264 * t207 - t225 * t249 + t318;
t314 = t249 * t206 - t207 * t248 + t317;
t300 = Icges(2,4) * t313;
t290 = rSges(2,1) * t313 - t310 * rSges(2,2);
t289 = t310 * rSges(2,1) + rSges(2,2) * t313;
t288 = Icges(2,1) * t313 - t351;
t287 = Icges(2,1) * t310 + t300;
t286 = -Icges(2,2) * t310 + t300;
t285 = Icges(2,2) * t313 + t351;
t278 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t277 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t276 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t259 = rSges(3,3) * t305 + (rSges(3,1) * t309 + rSges(3,2) * t312) * t303;
t256 = Icges(3,5) * t305 + (Icges(3,1) * t309 + Icges(3,4) * t312) * t303;
t254 = Icges(3,3) * t305 + (Icges(3,5) * t309 + Icges(3,6) * t312) * t303;
t253 = V_base(5) * rSges(2,3) - t289 * t299 + t333;
t252 = t290 * t299 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t250 = t289 * V_base(4) - t290 * V_base(5) + V_base(3);
t235 = qJD(5) * t257 + t264;
t234 = rSges(3,1) * t271 - rSges(3,2) * t270 + rSges(3,3) * t346;
t233 = t269 * rSges(3,1) - t268 * rSges(3,2) - rSges(3,3) * t344;
t232 = Icges(3,1) * t271 - Icges(3,4) * t270 + Icges(3,5) * t346;
t231 = Icges(3,1) * t269 - Icges(3,4) * t268 - Icges(3,5) * t344;
t228 = Icges(3,5) * t271 - Icges(3,6) * t270 + Icges(3,3) * t346;
t227 = Icges(3,5) * t269 - Icges(3,6) * t268 - Icges(3,3) * t344;
t224 = rSges(4,1) * t267 + rSges(4,2) * t266 - rSges(4,3) * t345;
t223 = Icges(4,1) * t267 + Icges(4,4) * t266 - Icges(4,5) * t345;
t222 = Icges(4,4) * t267 + Icges(4,2) * t266 - Icges(4,6) * t345;
t218 = rSges(5,1) * t258 - rSges(5,2) * t257 - rSges(5,3) * t345;
t217 = Icges(5,1) * t258 - Icges(5,4) * t257 - Icges(5,5) * t345;
t216 = Icges(5,4) * t258 - Icges(5,2) * t257 - Icges(5,6) * t345;
t215 = Icges(5,5) * t258 - Icges(5,6) * t257 - Icges(5,3) * t345;
t209 = qJD(5) * t242 + t249;
t208 = qJD(5) * t240 + t248;
t204 = rSges(4,1) * t247 + rSges(4,2) * t246 + rSges(4,3) * t270;
t203 = rSges(4,1) * t245 + rSges(4,2) * t244 + rSges(4,3) * t268;
t202 = Icges(4,1) * t247 + Icges(4,4) * t246 + Icges(4,5) * t270;
t201 = Icges(4,1) * t245 + Icges(4,4) * t244 + Icges(4,5) * t268;
t200 = Icges(4,4) * t247 + Icges(4,2) * t246 + Icges(4,6) * t270;
t199 = Icges(4,4) * t245 + Icges(4,2) * t244 + Icges(4,6) * t268;
t194 = rSges(5,1) * t243 - rSges(5,2) * t242 + rSges(5,3) * t270;
t193 = rSges(5,1) * t241 - rSges(5,2) * t240 + rSges(5,3) * t268;
t191 = Icges(5,1) * t243 - Icges(5,4) * t242 + Icges(5,5) * t270;
t190 = Icges(5,1) * t241 - Icges(5,4) * t240 + Icges(5,5) * t268;
t189 = Icges(5,4) * t243 - Icges(5,2) * t242 + Icges(5,6) * t270;
t188 = Icges(5,4) * t241 - Icges(5,2) * t240 + Icges(5,6) * t268;
t187 = Icges(5,5) * t243 - Icges(5,6) * t242 + Icges(5,3) * t270;
t186 = Icges(5,5) * t241 - Icges(5,6) * t240 + Icges(5,3) * t268;
t185 = rSges(6,1) * t239 + rSges(6,2) * t238 + rSges(6,3) * t257;
t173 = -t233 * t282 + t259 * t280 + t325;
t172 = t234 * t282 - t259 * t281 + t322;
t171 = t233 * t281 - t234 * t280 + t324;
t170 = rSges(6,1) * t214 + rSges(6,2) * t213 + rSges(6,3) * t242;
t168 = rSges(6,1) * t212 + rSges(6,2) * t211 + rSges(6,3) * t240;
t152 = t224 * t280 + (-t203 - t236) * t282 + t323;
t151 = t204 * t282 + (-t224 - t272) * t281 + t321;
t150 = t203 * t281 + (-t204 - t237) * t280 + t320;
t149 = -t193 * t264 + t218 * t248 + t319;
t148 = t194 * t264 - t218 * t249 + t318;
t147 = t193 * t249 - t194 * t248 + t317;
t146 = -t168 * t235 + t185 * t208 + t316;
t145 = t170 * t235 - t185 * t209 + t315;
t144 = t168 * t209 - t170 * t208 + t314;
t143 = qJD(6) * t242 + t208 * t336 - t235 * t338 + t316;
t142 = qJD(6) * t240 - t209 * t336 + t235 * t337 + t315;
t141 = qJD(6) * t257 - t208 * t337 + t209 * t338 + t314;
t1 = m(2) * (t250 ^ 2 + t252 ^ 2 + t253 ^ 2) / 0.2e1 + t248 * ((t187 * t268 - t189 * t240 + t191 * t241) * t249 + (t186 * t268 - t188 * t240 + t190 * t241) * t248 + (t215 * t268 - t216 * t240 + t217 * t241) * t264) / 0.2e1 + t249 * ((t187 * t270 - t189 * t242 + t191 * t243) * t249 + (t186 * t270 - t188 * t242 + t190 * t243) * t248 + (t215 * t270 - t216 * t242 + t217 * t243) * t264) / 0.2e1 + m(1) * (t276 ^ 2 + t277 ^ 2 + t278 ^ 2) / 0.2e1 + m(4) * (t150 ^ 2 + t151 ^ 2 + t152 ^ 2) / 0.2e1 + m(3) * (t171 ^ 2 + t172 ^ 2 + t173 ^ 2) / 0.2e1 + t264 * ((-t187 * t345 - t189 * t257 + t191 * t258) * t249 + (-t186 * t345 - t188 * t257 + t190 * t258) * t248 + (-t215 * t345 - t216 * t257 + t217 * t258) * t264) / 0.2e1 + m(5) * (t147 ^ 2 + t148 ^ 2 + t149 ^ 2) / 0.2e1 + m(6) * (t144 ^ 2 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + m(7) * (t141 ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + ((t211 * t366 + t212 * t365 + t240 * t367) * t235 + (t211 * t370 + t212 * t368 + t240 * t372) * t209 + (t371 * t211 + t369 * t212 + t373 * t240) * t208) * t208 / 0.2e1 + ((t213 * t366 + t214 * t365 + t242 * t367) * t235 + (t370 * t213 + t368 * t214 + t372 * t242) * t209 + (t371 * t213 + t369 * t214 + t242 * t373) * t208) * t209 / 0.2e1 + ((t366 * t238 + t365 * t239 + t367 * t257) * t235 + (t238 * t370 + t239 * t368 + t257 * t372) * t209 + (t371 * t238 + t369 * t239 + t257 * t373) * t208) * t235 / 0.2e1 + ((t222 * t244 + t223 * t245 - t254 * t344 + t269 * t256 + t268 * t362) * t282 + (t200 * t244 + t202 * t245 - t228 * t344 + t269 * t232 + t268 * t363) * t281 + (t199 * t244 + t201 * t245 - t227 * t344 + t269 * t231 + t364 * t268) * t280) * t280 / 0.2e1 + ((t222 * t246 + t223 * t247 + t254 * t346 + t256 * t271 + t270 * t362) * t282 + (t200 * t246 + t202 * t247 + t228 * t346 + t232 * t271 + t363 * t270) * t281 + (t199 * t246 + t201 * t247 + t227 * t346 + t231 * t271 + t270 * t364) * t280) * t281 / 0.2e1 + ((t227 * t280 + t228 * t281 + t254 * t282) * t305 + ((t230 * t312 + t232 * t309) * t281 + (t229 * t312 + t231 * t309) * t280 + (t255 * t312 + t256 * t309) * t282) * t303 + (-t198 * t345 + t200 * t266 + t202 * t267) * t281 + (-t197 * t345 + t199 * t266 + t201 * t267) * t280 + (-t221 * t345 + t222 * t266 + t223 * t267) * t282) * t282 / 0.2e1 + ((-t310 * t285 + t287 * t313 + Icges(1,4)) * V_base(5) + (-t310 * t286 + t288 * t313 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t285 * t313 + t310 * t287 + Icges(1,2)) * V_base(5) + (t286 * t313 + t310 * t288 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t310 + Icges(2,6) * t313) * V_base(5) + (Icges(2,5) * t313 - Icges(2,6) * t310) * V_base(4) + Icges(2,3) * t299 / 0.2e1) * t299;
T  = t1;
