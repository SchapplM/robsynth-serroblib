% Calculate kinetic energy for
% S6RRRPRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 17:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRP6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP6_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP6_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPRP6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP6_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP6_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP6_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRP6_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:53:47
% EndTime: 2019-03-09 16:53:51
% DurationCPUTime: 4.42s
% Computational Cost: add. (3512->386), mult. (5968->541), div. (0->0), fcn. (7093->12), ass. (0->175)
t382 = Icges(6,1) + Icges(7,1);
t381 = Icges(6,4) + Icges(7,4);
t380 = Icges(6,5) + Icges(7,5);
t379 = Icges(6,2) + Icges(7,2);
t378 = Icges(6,6) + Icges(7,6);
t377 = Icges(4,3) + Icges(5,3);
t376 = Icges(6,3) + Icges(7,3);
t375 = rSges(7,3) + qJ(6);
t303 = cos(pkin(6));
t309 = sin(qJ(1));
t312 = cos(qJ(2));
t340 = t309 * t312;
t308 = sin(qJ(2));
t313 = cos(qJ(1));
t341 = t308 * t313;
t269 = t303 * t341 + t340;
t334 = qJ(3) + pkin(11);
t298 = sin(t334);
t327 = cos(t334);
t302 = sin(pkin(6));
t344 = t302 * t313;
t241 = t269 * t327 - t298 * t344;
t339 = t312 * t313;
t342 = t308 * t309;
t268 = -t303 * t339 + t342;
t306 = sin(qJ(5));
t310 = cos(qJ(5));
t211 = -t241 * t306 + t268 * t310;
t350 = t268 * t306;
t212 = t241 * t310 + t350;
t326 = t302 * t327;
t240 = t269 * t298 + t313 * t326;
t374 = t378 * t211 + t380 * t212 + t376 * t240;
t271 = -t303 * t342 + t339;
t347 = t302 * t309;
t243 = t271 * t327 + t298 * t347;
t270 = t303 * t340 + t341;
t213 = -t243 * t306 + t270 * t310;
t349 = t270 * t306;
t214 = t243 * t310 + t349;
t242 = t271 * t298 - t309 * t326;
t373 = t378 * t213 + t380 * t214 + t376 * t242;
t372 = t379 * t211 + t381 * t212 + t378 * t240;
t371 = t379 * t213 + t381 * t214 + t378 * t242;
t370 = t381 * t211 + t382 * t212 + t380 * t240;
t369 = t381 * t213 + t382 * t214 + t380 * t242;
t258 = t303 * t298 + t308 * t326;
t345 = t302 * t312;
t238 = -t258 * t306 - t310 * t345;
t330 = t306 * t345;
t239 = t258 * t310 - t330;
t348 = t302 * t308;
t257 = t298 * t348 - t303 * t327;
t368 = t378 * t238 + t380 * t239 + t376 * t257;
t367 = t379 * t238 + t381 * t239 + t378 * t257;
t366 = t381 * t238 + t382 * t239 + t380 * t257;
t307 = sin(qJ(3));
t311 = cos(qJ(3));
t246 = -t269 * t307 - t311 * t344;
t328 = t307 * t344;
t247 = t269 * t311 - t328;
t365 = Icges(4,5) * t247 + Icges(5,5) * t241 + Icges(4,6) * t246 - Icges(5,6) * t240 + t377 * t268;
t346 = t302 * t311;
t248 = -t271 * t307 + t309 * t346;
t329 = t307 * t347;
t249 = t271 * t311 + t329;
t364 = Icges(4,5) * t249 + Icges(5,5) * t243 + Icges(4,6) * t248 - Icges(5,6) * t242 + t377 * t270;
t266 = t303 * t311 - t307 * t348;
t343 = t303 * t307;
t267 = t308 * t346 + t343;
t363 = Icges(4,5) * t267 + Icges(5,5) * t258 + Icges(4,6) * t266 - Icges(5,6) * t257 - t377 * t345;
t356 = pkin(8) * t303;
t355 = pkin(3) * t311;
t354 = pkin(5) * t310;
t351 = Icges(2,4) * t309;
t338 = rSges(7,1) * t212 + rSges(7,2) * t211 + pkin(5) * t350 + t240 * t375 + t241 * t354;
t337 = rSges(7,1) * t214 + rSges(7,2) * t213 + pkin(5) * t349 + t242 * t375 + t243 * t354;
t336 = rSges(7,1) * t239 + rSges(7,2) * t238 - pkin(5) * t330 + t257 * t375 + t258 * t354;
t335 = qJD(2) * t302;
t333 = V_base(5) * pkin(7) + V_base(1);
t281 = t309 * t335 + V_base(4);
t299 = V_base(6) + qJD(1);
t245 = qJD(3) * t270 + t281;
t282 = qJD(2) * t303 + t299;
t280 = -t313 * t335 + V_base(5);
t274 = t309 * pkin(1) - pkin(8) * t344;
t325 = -t274 * t299 + V_base(5) * t356 + t333;
t275 = pkin(1) * t313 + pkin(8) * t347;
t324 = V_base(4) * t274 - t275 * V_base(5) + V_base(3);
t244 = qJD(3) * t268 + t280;
t264 = -qJD(3) * t345 + t282;
t323 = t299 * t275 + V_base(2) + (-pkin(7) - t356) * V_base(4);
t236 = t269 * pkin(2) + t268 * pkin(9);
t273 = (pkin(2) * t308 - pkin(9) * t312) * t302;
t322 = -t236 * t282 + t280 * t273 + t325;
t237 = pkin(2) * t271 + pkin(9) * t270;
t321 = t281 * t236 - t237 * t280 + t324;
t225 = pkin(3) * t343 + (-qJ(4) * t312 + t308 * t355) * t302;
t320 = qJD(4) * t270 + t244 * t225 + t322;
t319 = t282 * t237 - t273 * t281 + t323;
t196 = pkin(3) * t329 + qJ(4) * t270 + t271 * t355;
t318 = qJD(4) * t268 + t264 * t196 + t319;
t195 = -pkin(3) * t328 + qJ(4) * t268 + t269 * t355;
t317 = -qJD(4) * t345 + t245 * t195 + t321;
t207 = pkin(4) * t241 + pkin(10) * t240;
t221 = pkin(4) * t258 + pkin(10) * t257;
t316 = t244 * t221 + (-t195 - t207) * t264 + t320;
t208 = pkin(4) * t243 + pkin(10) * t242;
t315 = t264 * t208 + (-t221 - t225) * t245 + t318;
t314 = t245 * t207 + (-t196 - t208) * t244 + t317;
t300 = Icges(2,4) * t313;
t290 = rSges(2,1) * t313 - t309 * rSges(2,2);
t289 = t309 * rSges(2,1) + rSges(2,2) * t313;
t288 = Icges(2,1) * t313 - t351;
t287 = Icges(2,1) * t309 + t300;
t286 = -Icges(2,2) * t309 + t300;
t285 = Icges(2,2) * t313 + t351;
t278 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t277 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t276 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t259 = rSges(3,3) * t303 + (rSges(3,1) * t308 + rSges(3,2) * t312) * t302;
t256 = Icges(3,5) * t303 + (Icges(3,1) * t308 + Icges(3,4) * t312) * t302;
t255 = Icges(3,6) * t303 + (Icges(3,4) * t308 + Icges(3,2) * t312) * t302;
t254 = Icges(3,3) * t303 + (Icges(3,5) * t308 + Icges(3,6) * t312) * t302;
t253 = V_base(5) * rSges(2,3) - t289 * t299 + t333;
t252 = t290 * t299 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t250 = t289 * V_base(4) - t290 * V_base(5) + V_base(3);
t235 = qJD(5) * t257 + t264;
t234 = rSges(3,1) * t271 - rSges(3,2) * t270 + rSges(3,3) * t347;
t233 = t269 * rSges(3,1) - t268 * rSges(3,2) - rSges(3,3) * t344;
t232 = Icges(3,1) * t271 - Icges(3,4) * t270 + Icges(3,5) * t347;
t231 = Icges(3,1) * t269 - Icges(3,4) * t268 - Icges(3,5) * t344;
t230 = Icges(3,4) * t271 - Icges(3,2) * t270 + Icges(3,6) * t347;
t229 = Icges(3,4) * t269 - Icges(3,2) * t268 - Icges(3,6) * t344;
t228 = Icges(3,5) * t271 - Icges(3,6) * t270 + Icges(3,3) * t347;
t227 = Icges(3,5) * t269 - Icges(3,6) * t268 - Icges(3,3) * t344;
t226 = rSges(4,1) * t267 + rSges(4,2) * t266 - rSges(4,3) * t345;
t224 = Icges(4,1) * t267 + Icges(4,4) * t266 - Icges(4,5) * t345;
t223 = Icges(4,4) * t267 + Icges(4,2) * t266 - Icges(4,6) * t345;
t218 = rSges(5,1) * t258 - rSges(5,2) * t257 - rSges(5,3) * t345;
t217 = Icges(5,1) * t258 - Icges(5,4) * t257 - Icges(5,5) * t345;
t216 = Icges(5,4) * t258 - Icges(5,2) * t257 - Icges(5,6) * t345;
t210 = qJD(5) * t242 + t245;
t209 = qJD(5) * t240 + t244;
t204 = rSges(4,1) * t249 + rSges(4,2) * t248 + rSges(4,3) * t270;
t203 = rSges(4,1) * t247 + rSges(4,2) * t246 + rSges(4,3) * t268;
t202 = Icges(4,1) * t249 + Icges(4,4) * t248 + Icges(4,5) * t270;
t201 = Icges(4,1) * t247 + Icges(4,4) * t246 + Icges(4,5) * t268;
t200 = Icges(4,4) * t249 + Icges(4,2) * t248 + Icges(4,6) * t270;
t199 = Icges(4,4) * t247 + Icges(4,2) * t246 + Icges(4,6) * t268;
t194 = rSges(5,1) * t243 - rSges(5,2) * t242 + rSges(5,3) * t270;
t193 = rSges(5,1) * t241 - rSges(5,2) * t240 + rSges(5,3) * t268;
t191 = Icges(5,1) * t243 - Icges(5,4) * t242 + Icges(5,5) * t270;
t190 = Icges(5,1) * t241 - Icges(5,4) * t240 + Icges(5,5) * t268;
t189 = Icges(5,4) * t243 - Icges(5,2) * t242 + Icges(5,6) * t270;
t188 = Icges(5,4) * t241 - Icges(5,2) * t240 + Icges(5,6) * t268;
t185 = rSges(6,1) * t239 + rSges(6,2) * t238 + rSges(6,3) * t257;
t173 = -t233 * t282 + t259 * t280 + t325;
t172 = t234 * t282 - t259 * t281 + t323;
t171 = t233 * t281 - t234 * t280 + t324;
t170 = rSges(6,1) * t214 + rSges(6,2) * t213 + rSges(6,3) * t242;
t168 = rSges(6,1) * t212 + rSges(6,2) * t211 + rSges(6,3) * t240;
t152 = -t203 * t264 + t226 * t244 + t322;
t151 = t204 * t264 - t226 * t245 + t319;
t150 = t203 * t245 - t204 * t244 + t321;
t149 = t218 * t244 + (-t193 - t195) * t264 + t320;
t148 = t194 * t264 + (-t218 - t225) * t245 + t318;
t147 = t193 * t245 + (-t194 - t196) * t244 + t317;
t146 = -t168 * t235 + t185 * t209 + t316;
t145 = t170 * t235 - t185 * t210 + t315;
t144 = t168 * t210 - t170 * t209 + t314;
t143 = qJD(6) * t242 + t209 * t336 - t235 * t338 + t316;
t142 = qJD(6) * t240 - t210 * t336 + t235 * t337 + t315;
t141 = qJD(6) * t257 - t209 * t337 + t210 * t338 + t314;
t1 = t280 * ((-t228 * t344 - t268 * t230 + t269 * t232) * t281 + (-t227 * t344 - t268 * t229 + t269 * t231) * t280 + (-t254 * t344 - t268 * t255 + t269 * t256) * t282) / 0.2e1 + t281 * ((t228 * t347 - t230 * t270 + t232 * t271) * t281 + (t227 * t347 - t229 * t270 + t231 * t271) * t280 + (t254 * t347 - t255 * t270 + t256 * t271) * t282) / 0.2e1 + m(4) * (t150 ^ 2 + t151 ^ 2 + t152 ^ 2) / 0.2e1 + m(5) * (t147 ^ 2 + t148 ^ 2 + t149 ^ 2) / 0.2e1 + m(6) * (t144 ^ 2 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + m(7) * (t141 ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + t282 * ((t227 * t280 + t228 * t281 + t254 * t282) * t303 + ((t230 * t312 + t232 * t308) * t281 + (t229 * t312 + t231 * t308) * t280 + (t255 * t312 + t256 * t308) * t282) * t302) / 0.2e1 + m(3) * (t171 ^ 2 + t172 ^ 2 + t173 ^ 2) / 0.2e1 + m(2) * (t250 ^ 2 + t252 ^ 2 + t253 ^ 2) / 0.2e1 + m(1) * (t276 ^ 2 + t277 ^ 2 + t278 ^ 2) / 0.2e1 + ((t211 * t367 + t212 * t366 + t240 * t368) * t235 + (t211 * t371 + t212 * t369 + t240 * t373) * t210 + (t372 * t211 + t370 * t212 + t374 * t240) * t209) * t209 / 0.2e1 + ((t213 * t367 + t214 * t366 + t242 * t368) * t235 + (t371 * t213 + t369 * t214 + t373 * t242) * t210 + (t372 * t213 + t370 * t214 + t242 * t374) * t209) * t210 / 0.2e1 + ((t367 * t238 + t366 * t239 + t368 * t257) * t235 + (t238 * t371 + t239 * t369 + t257 * t373) * t210 + (t372 * t238 + t370 * t239 + t257 * t374) * t209) * t235 / 0.2e1 + ((-t216 * t240 + t217 * t241 + t223 * t246 + t224 * t247 + t268 * t363) * t264 + (-t189 * t240 + t191 * t241 + t200 * t246 + t202 * t247 + t268 * t364) * t245 + (-t188 * t240 + t190 * t241 + t199 * t246 + t201 * t247 + t365 * t268) * t244) * t244 / 0.2e1 + ((-t216 * t242 + t217 * t243 + t223 * t248 + t224 * t249 + t270 * t363) * t264 + (-t189 * t242 + t191 * t243 + t200 * t248 + t202 * t249 + t364 * t270) * t245 + (-t188 * t242 + t190 * t243 + t199 * t248 + t201 * t249 + t270 * t365) * t244) * t245 / 0.2e1 + ((-t216 * t257 + t217 * t258 + t223 * t266 + t224 * t267 - t363 * t345) * t264 + (-t189 * t257 + t191 * t258 + t200 * t266 + t202 * t267 - t345 * t364) * t245 + (-t188 * t257 + t190 * t258 + t199 * t266 + t201 * t267 - t345 * t365) * t244) * t264 / 0.2e1 + ((-t309 * t285 + t287 * t313 + Icges(1,4)) * V_base(5) + (-t309 * t286 + t288 * t313 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t285 * t313 + t309 * t287 + Icges(1,2)) * V_base(5) + (t286 * t313 + t309 * t288 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t309 + Icges(2,6) * t313) * V_base(5) + (Icges(2,5) * t313 - Icges(2,6) * t309) * V_base(4) + Icges(2,3) * t299 / 0.2e1) * t299;
T  = t1;
