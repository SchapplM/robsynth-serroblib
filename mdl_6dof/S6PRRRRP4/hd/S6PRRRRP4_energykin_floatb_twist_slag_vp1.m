% Calculate kinetic energy for
% S6PRRRRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRRP4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP4_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRRRP4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP4_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRP4_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRRP4_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:13:01
% EndTime: 2019-03-09 00:13:06
% DurationCPUTime: 4.42s
% Computational Cost: add. (3320->386), mult. (6925->560), div. (0->0), fcn. (8407->12), ass. (0->174)
t385 = Icges(6,1) + Icges(7,1);
t384 = -Icges(6,4) + Icges(7,5);
t383 = Icges(7,4) + Icges(6,5);
t382 = Icges(6,2) + Icges(7,3);
t381 = Icges(7,2) + Icges(6,3);
t380 = -Icges(6,6) + Icges(7,6);
t379 = rSges(7,1) + pkin(5);
t378 = rSges(7,3) + qJ(6);
t316 = sin(pkin(11));
t318 = cos(pkin(11));
t324 = cos(qJ(2));
t319 = cos(pkin(6));
t322 = sin(qJ(2));
t352 = t319 * t322;
t281 = t316 * t324 + t318 * t352;
t317 = sin(pkin(6));
t321 = sin(qJ(3));
t354 = t317 * t321;
t363 = cos(qJ(3));
t263 = t281 * t363 - t318 * t354;
t351 = t319 * t324;
t280 = t316 * t322 - t318 * t351;
t350 = qJ(4) + qJ(5);
t314 = sin(t350);
t338 = cos(t350);
t228 = t263 * t314 - t280 * t338;
t229 = t263 * t338 + t280 * t314;
t339 = t317 * t363;
t262 = t281 * t321 + t318 * t339;
t376 = t228 * t382 + t229 * t384 + t262 * t380;
t283 = -t316 * t352 + t318 * t324;
t265 = t283 * t363 + t316 * t354;
t282 = t316 * t351 + t318 * t322;
t230 = t265 * t314 - t282 * t338;
t231 = t265 * t338 + t282 * t314;
t264 = t283 * t321 - t316 * t339;
t375 = t230 * t382 + t231 * t384 + t264 * t380;
t374 = t228 * t380 + t229 * t383 + t262 * t381;
t373 = t230 * t380 + t231 * t383 + t264 * t381;
t372 = t228 * t384 + t229 * t385 + t262 * t383;
t371 = t230 * t384 + t231 * t385 + t264 * t383;
t288 = t319 * t321 + t322 * t339;
t353 = t317 * t324;
t255 = t288 * t314 + t338 * t353;
t256 = t288 * t338 - t314 * t353;
t287 = -t319 * t363 + t322 * t354;
t370 = t255 * t382 + t256 * t384 + t287 * t380;
t369 = t255 * t380 + t256 * t383 + t287 * t381;
t368 = t255 * t384 + t256 * t385 + t287 * t383;
t362 = pkin(7) * t319;
t323 = cos(qJ(4));
t361 = pkin(4) * t323;
t359 = Icges(2,4) * t316;
t320 = sin(qJ(4));
t358 = t280 * t320;
t357 = t282 * t320;
t356 = t316 * t317;
t355 = t317 * t318;
t349 = rSges(7,2) * t262 + t228 * t378 + t229 * t379;
t348 = rSges(7,2) * t264 + t230 * t378 + t231 * t379;
t347 = rSges(7,2) * t287 + t255 * t378 + t256 * t379;
t346 = qJD(2) * t317;
t345 = V_base(5) * qJ(1) + V_base(1);
t341 = qJD(1) + V_base(3);
t340 = t320 * t353;
t296 = t316 * t346 + V_base(4);
t307 = qJD(2) * t319 + V_base(6);
t261 = qJD(3) * t282 + t296;
t227 = qJD(4) * t264 + t261;
t295 = -t318 * t346 + V_base(5);
t260 = qJD(3) * t280 + t295;
t284 = -qJD(3) * t353 + t307;
t290 = pkin(1) * t316 - pkin(7) * t355;
t337 = -t290 * V_base(6) + t362 * V_base(5) + t345;
t291 = pkin(1) * t318 + pkin(7) * t356;
t336 = t290 * V_base(4) - t291 * V_base(5) + t341;
t226 = qJD(4) * t262 + t260;
t254 = qJD(4) * t287 + t284;
t335 = V_base(6) * t291 + V_base(2) + (-qJ(1) - t362) * V_base(4);
t251 = pkin(2) * t281 + pkin(8) * t280;
t289 = (pkin(2) * t322 - pkin(8) * t324) * t317;
t334 = -t251 * t307 + t289 * t295 + t337;
t252 = pkin(2) * t283 + pkin(8) * t282;
t333 = t251 * t296 - t252 * t295 + t336;
t332 = t252 * t307 - t289 * t296 + t335;
t224 = pkin(3) * t263 + pkin(9) * t262;
t253 = pkin(3) * t288 + pkin(9) * t287;
t331 = -t224 * t284 + t253 * t260 + t334;
t225 = pkin(3) * t265 + pkin(9) * t264;
t330 = t224 * t261 - t225 * t260 + t333;
t329 = t225 * t284 - t253 * t261 + t332;
t165 = pkin(4) * t358 + pkin(10) * t262 + t263 * t361;
t208 = -pkin(4) * t340 + pkin(10) * t287 + t288 * t361;
t328 = -t165 * t254 + t208 * t226 + t331;
t166 = pkin(4) * t357 + pkin(10) * t264 + t265 * t361;
t327 = t165 * t227 - t166 * t226 + t330;
t326 = t166 * t254 - t208 * t227 + t329;
t313 = Icges(2,4) * t318;
t304 = rSges(2,1) * t318 - rSges(2,2) * t316;
t303 = rSges(2,1) * t316 + rSges(2,2) * t318;
t302 = Icges(2,1) * t318 - t359;
t301 = Icges(2,1) * t316 + t313;
t300 = -Icges(2,2) * t316 + t313;
t299 = Icges(2,2) * t318 + t359;
t294 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t293 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t292 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t276 = rSges(3,3) * t319 + (rSges(3,1) * t322 + rSges(3,2) * t324) * t317;
t275 = Icges(3,5) * t319 + (Icges(3,1) * t322 + Icges(3,4) * t324) * t317;
t274 = Icges(3,6) * t319 + (Icges(3,4) * t322 + Icges(3,2) * t324) * t317;
t273 = Icges(3,3) * t319 + (Icges(3,5) * t322 + Icges(3,6) * t324) * t317;
t270 = V_base(5) * rSges(2,3) - t303 * V_base(6) + t345;
t269 = t304 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t267 = t288 * t323 - t340;
t266 = -t288 * t320 - t323 * t353;
t259 = t303 * V_base(4) - t304 * V_base(5) + t341;
t250 = rSges(4,1) * t288 - rSges(4,2) * t287 - rSges(4,3) * t353;
t249 = Icges(4,1) * t288 - Icges(4,4) * t287 - Icges(4,5) * t353;
t248 = Icges(4,4) * t288 - Icges(4,2) * t287 - Icges(4,6) * t353;
t247 = Icges(4,5) * t288 - Icges(4,6) * t287 - Icges(4,3) * t353;
t246 = rSges(3,1) * t283 - rSges(3,2) * t282 + rSges(3,3) * t356;
t245 = rSges(3,1) * t281 - rSges(3,2) * t280 - rSges(3,3) * t355;
t244 = Icges(3,1) * t283 - Icges(3,4) * t282 + Icges(3,5) * t356;
t243 = Icges(3,1) * t281 - Icges(3,4) * t280 - Icges(3,5) * t355;
t242 = Icges(3,4) * t283 - Icges(3,2) * t282 + Icges(3,6) * t356;
t241 = Icges(3,4) * t281 - Icges(3,2) * t280 - Icges(3,6) * t355;
t240 = Icges(3,5) * t283 - Icges(3,6) * t282 + Icges(3,3) * t356;
t239 = Icges(3,5) * t281 - Icges(3,6) * t280 - Icges(3,3) * t355;
t236 = qJD(5) * t287 + t254;
t235 = t265 * t323 + t357;
t234 = -t265 * t320 + t282 * t323;
t233 = t263 * t323 + t358;
t232 = -t263 * t320 + t280 * t323;
t221 = rSges(5,1) * t267 + rSges(5,2) * t266 + rSges(5,3) * t287;
t219 = Icges(5,1) * t267 + Icges(5,4) * t266 + Icges(5,5) * t287;
t218 = Icges(5,4) * t267 + Icges(5,2) * t266 + Icges(5,6) * t287;
t217 = Icges(5,5) * t267 + Icges(5,6) * t266 + Icges(5,3) * t287;
t216 = rSges(4,1) * t265 - rSges(4,2) * t264 + rSges(4,3) * t282;
t215 = rSges(4,1) * t263 - rSges(4,2) * t262 + rSges(4,3) * t280;
t214 = Icges(4,1) * t265 - Icges(4,4) * t264 + Icges(4,5) * t282;
t213 = Icges(4,1) * t263 - Icges(4,4) * t262 + Icges(4,5) * t280;
t212 = Icges(4,4) * t265 - Icges(4,2) * t264 + Icges(4,6) * t282;
t211 = Icges(4,4) * t263 - Icges(4,2) * t262 + Icges(4,6) * t280;
t210 = Icges(4,5) * t265 - Icges(4,6) * t264 + Icges(4,3) * t282;
t209 = Icges(4,5) * t263 - Icges(4,6) * t262 + Icges(4,3) * t280;
t207 = rSges(6,1) * t256 - rSges(6,2) * t255 + rSges(6,3) * t287;
t199 = qJD(5) * t264 + t227;
t198 = qJD(5) * t262 + t226;
t196 = -t245 * t307 + t276 * t295 + t337;
t195 = t246 * t307 - t276 * t296 + t335;
t192 = rSges(5,1) * t235 + rSges(5,2) * t234 + rSges(5,3) * t264;
t191 = rSges(5,1) * t233 + rSges(5,2) * t232 + rSges(5,3) * t262;
t190 = Icges(5,1) * t235 + Icges(5,4) * t234 + Icges(5,5) * t264;
t189 = Icges(5,1) * t233 + Icges(5,4) * t232 + Icges(5,5) * t262;
t188 = Icges(5,4) * t235 + Icges(5,2) * t234 + Icges(5,6) * t264;
t187 = Icges(5,4) * t233 + Icges(5,2) * t232 + Icges(5,6) * t262;
t186 = Icges(5,5) * t235 + Icges(5,6) * t234 + Icges(5,3) * t264;
t185 = Icges(5,5) * t233 + Icges(5,6) * t232 + Icges(5,3) * t262;
t183 = t245 * t296 - t246 * t295 + t336;
t182 = rSges(6,1) * t231 - rSges(6,2) * t230 + rSges(6,3) * t264;
t180 = rSges(6,1) * t229 - rSges(6,2) * t228 + rSges(6,3) * t262;
t162 = -t215 * t284 + t250 * t260 + t334;
t161 = t216 * t284 - t250 * t261 + t332;
t160 = t215 * t261 - t216 * t260 + t333;
t159 = -t191 * t254 + t221 * t226 + t331;
t158 = t192 * t254 - t221 * t227 + t329;
t157 = t191 * t227 - t192 * t226 + t330;
t156 = -t180 * t236 + t198 * t207 + t328;
t155 = t182 * t236 - t199 * t207 + t326;
t154 = t180 * t199 - t182 * t198 + t327;
t153 = qJD(6) * t230 + t198 * t347 - t236 * t349 + t328;
t152 = qJD(6) * t228 - t199 * t347 + t236 * t348 + t326;
t151 = qJD(6) * t255 - t198 * t348 + t199 * t349 + t327;
t1 = t295 * ((-t240 * t355 - t242 * t280 + t244 * t281) * t296 + (-t239 * t355 - t241 * t280 + t243 * t281) * t295 + (-t273 * t355 - t274 * t280 + t275 * t281) * t307) / 0.2e1 + t296 * ((t240 * t356 - t242 * t282 + t244 * t283) * t296 + (t239 * t356 - t241 * t282 + t243 * t283) * t295 + (t273 * t356 - t274 * t282 + t275 * t283) * t307) / 0.2e1 + t284 * ((-t210 * t353 - t212 * t287 + t214 * t288) * t261 + (-t209 * t353 - t211 * t287 + t213 * t288) * t260 + (-t247 * t353 - t248 * t287 + t249 * t288) * t284) / 0.2e1 + m(1) * (t292 ^ 2 + t293 ^ 2 + t294 ^ 2) / 0.2e1 + m(7) * (t151 ^ 2 + t152 ^ 2 + t153 ^ 2) / 0.2e1 + t260 * ((t210 * t280 - t212 * t262 + t214 * t263) * t261 + (t209 * t280 - t211 * t262 + t213 * t263) * t260 + (t247 * t280 - t248 * t262 + t249 * t263) * t284) / 0.2e1 + t254 * ((t186 * t287 + t188 * t266 + t190 * t267) * t227 + (t185 * t287 + t187 * t266 + t189 * t267) * t226 + (t217 * t287 + t218 * t266 + t219 * t267) * t254) / 0.2e1 + t261 * ((t210 * t282 - t212 * t264 + t214 * t265) * t261 + (t209 * t282 - t211 * t264 + t213 * t265) * t260 + (t247 * t282 - t248 * t264 + t249 * t265) * t284) / 0.2e1 + m(2) * (t259 ^ 2 + t269 ^ 2 + t270 ^ 2) / 0.2e1 + t227 * ((t186 * t264 + t188 * t234 + t190 * t235) * t227 + (t185 * t264 + t187 * t234 + t189 * t235) * t226 + (t217 * t264 + t218 * t234 + t219 * t235) * t254) / 0.2e1 + t226 * ((t186 * t262 + t188 * t232 + t190 * t233) * t227 + (t185 * t262 + t187 * t232 + t189 * t233) * t226 + (t217 * t262 + t218 * t232 + t219 * t233) * t254) / 0.2e1 + t307 * ((t239 * t295 + t240 * t296 + t273 * t307) * t319 + ((t242 * t324 + t244 * t322) * t296 + (t241 * t324 + t243 * t322) * t295 + (t274 * t324 + t275 * t322) * t307) * t317) / 0.2e1 + m(3) * (t183 ^ 2 + t195 ^ 2 + t196 ^ 2) / 0.2e1 + m(4) * (t160 ^ 2 + t161 ^ 2 + t162 ^ 2) / 0.2e1 + m(5) * (t157 ^ 2 + t158 ^ 2 + t159 ^ 2) / 0.2e1 + m(6) * (t154 ^ 2 + t155 ^ 2 + t156 ^ 2) / 0.2e1 + ((t228 * t370 + t229 * t368 + t262 * t369) * t236 + (t228 * t375 + t229 * t371 + t262 * t373) * t199 + (t228 * t376 + t229 * t372 + t262 * t374) * t198) * t198 / 0.2e1 + ((t230 * t370 + t231 * t368 + t264 * t369) * t236 + (t230 * t375 + t231 * t371 + t264 * t373) * t199 + (t230 * t376 + t231 * t372 + t264 * t374) * t198) * t199 / 0.2e1 + ((t255 * t370 + t256 * t368 + t287 * t369) * t236 + (t255 * t375 + t256 * t371 + t287 * t373) * t199 + (t255 * t376 + t256 * t372 + t287 * t374) * t198) * t236 / 0.2e1 + ((-t299 * t316 + t301 * t318 + Icges(1,4)) * V_base(5) + (-t300 * t316 + t302 * t318 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t299 * t318 + t301 * t316 + Icges(1,2)) * V_base(5) + (t300 * t318 + t302 * t316 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t318 - Icges(2,6) * t316 + Icges(1,5)) * V_base(4) + (Icges(2,5) * t316 + Icges(2,6) * t318 + Icges(1,6)) * V_base(5) + (Icges(2,3) / 0.2e1 + Icges(1,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
