% Calculate kinetic energy for
% S6PRPRRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-03-08 20:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRRP4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP4_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRPRRP4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP4_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRP4_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRP4_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:09:09
% EndTime: 2019-03-08 20:09:13
% DurationCPUTime: 4.35s
% Computational Cost: add. (3344->386), mult. (5819->540), div. (0->0), fcn. (6923->12), ass. (0->174)
t383 = Icges(6,1) + Icges(7,1);
t382 = -Icges(6,4) + Icges(7,5);
t381 = Icges(7,4) + Icges(6,5);
t380 = Icges(6,2) + Icges(7,3);
t379 = Icges(7,2) + Icges(6,3);
t378 = -Icges(6,6) + Icges(7,6);
t377 = rSges(7,1) + pkin(5);
t376 = rSges(7,3) + qJ(6);
t309 = sin(pkin(10));
t312 = cos(pkin(10));
t317 = cos(qJ(2));
t313 = cos(pkin(6));
t316 = sin(qJ(2));
t346 = t313 * t316;
t275 = t309 * t317 + t312 * t346;
t339 = pkin(11) + qJ(4);
t305 = sin(t339);
t331 = cos(t339);
t310 = sin(pkin(6));
t349 = t310 * t312;
t245 = t275 * t331 - t305 * t349;
t345 = t313 * t317;
t274 = t309 * t316 - t312 * t345;
t315 = sin(qJ(5));
t355 = cos(qJ(5));
t216 = t245 * t315 - t274 * t355;
t217 = t245 * t355 + t274 * t315;
t330 = t310 * t331;
t244 = t275 * t305 + t312 * t330;
t374 = t380 * t216 + t382 * t217 + t378 * t244;
t277 = -t309 * t346 + t312 * t317;
t350 = t309 * t310;
t247 = t277 * t331 + t305 * t350;
t276 = t309 * t345 + t312 * t316;
t218 = t247 * t315 - t276 * t355;
t219 = t247 * t355 + t276 * t315;
t246 = t277 * t305 - t309 * t330;
t373 = t380 * t218 + t382 * t219 + t378 * t246;
t372 = t378 * t216 + t381 * t217 + t379 * t244;
t371 = t378 * t218 + t381 * t219 + t379 * t246;
t370 = t382 * t216 + t383 * t217 + t381 * t244;
t369 = t382 * t218 + t383 * t219 + t381 * t246;
t264 = t313 * t305 + t316 * t330;
t347 = t310 * t317;
t248 = t264 * t315 + t347 * t355;
t249 = t264 * t355 - t315 * t347;
t348 = t310 * t316;
t263 = t305 * t348 - t313 * t331;
t368 = t380 * t248 + t382 * t249 + t378 * t263;
t367 = t378 * t248 + t381 * t249 + t379 * t263;
t366 = t382 * t248 + t383 * t249 + t381 * t263;
t308 = sin(pkin(11));
t311 = cos(pkin(11));
t251 = -t275 * t308 - t311 * t349;
t332 = t308 * t349;
t252 = t275 * t311 - t332;
t202 = Icges(4,5) * t252 + Icges(4,6) * t251 + Icges(4,3) * t274;
t232 = Icges(3,4) * t275 - Icges(3,2) * t274 - Icges(3,6) * t349;
t365 = t202 - t232;
t253 = -t277 * t308 + t311 * t350;
t333 = t308 * t350;
t254 = t277 * t311 + t333;
t203 = Icges(4,5) * t254 + Icges(4,6) * t253 + Icges(4,3) * t276;
t233 = Icges(3,4) * t277 - Icges(3,2) * t276 + Icges(3,6) * t350;
t364 = t203 - t233;
t272 = -t308 * t348 + t311 * t313;
t351 = t308 * t313;
t273 = t311 * t348 + t351;
t227 = Icges(4,5) * t273 + Icges(4,6) * t272 - Icges(4,3) * t347;
t261 = Icges(3,6) * t313 + (Icges(3,4) * t316 + Icges(3,2) * t317) * t310;
t363 = t227 - t261;
t354 = pkin(7) * t313;
t353 = pkin(3) * t311;
t352 = Icges(2,4) * t309;
t343 = rSges(7,2) * t244 + t376 * t216 + t377 * t217;
t342 = rSges(7,2) * t246 + t376 * t218 + t377 * t219;
t341 = rSges(7,2) * t263 + t376 * t248 + t377 * t249;
t340 = qJD(2) * t310;
t338 = V_base(5) * qJ(1) + V_base(1);
t334 = qJD(1) + V_base(3);
t289 = t309 * t340 + V_base(4);
t300 = qJD(2) * t313 + V_base(6);
t256 = qJD(4) * t276 + t289;
t288 = -t312 * t340 + V_base(5);
t255 = qJD(4) * t274 + t288;
t278 = -qJD(4) * t347 + t300;
t282 = pkin(1) * t309 - pkin(7) * t349;
t329 = -t282 * V_base(6) + V_base(5) * t354 + t338;
t283 = pkin(1) * t312 + pkin(7) * t350;
t328 = V_base(4) * t282 - V_base(5) * t283 + t334;
t327 = V_base(6) * t283 + V_base(2) + (-qJ(1) - t354) * V_base(4);
t281 = (pkin(2) * t316 - qJ(3) * t317) * t310;
t326 = qJD(3) * t276 + t288 * t281 + t329;
t243 = pkin(2) * t277 + qJ(3) * t276;
t325 = qJD(3) * t274 + t300 * t243 + t327;
t242 = pkin(2) * t275 + qJ(3) * t274;
t324 = -qJD(3) * t347 + t289 * t242 + t328;
t199 = -pkin(3) * t332 + pkin(8) * t274 + t275 * t353;
t238 = pkin(3) * t351 + (-pkin(8) * t317 + t316 * t353) * t310;
t323 = t288 * t238 + (-t199 - t242) * t300 + t326;
t200 = pkin(3) * t333 + pkin(8) * t276 + t277 * t353;
t322 = t300 * t200 + (-t238 - t281) * t289 + t325;
t321 = t289 * t199 + (-t200 - t243) * t288 + t324;
t211 = pkin(4) * t245 + pkin(9) * t244;
t237 = pkin(4) * t264 + pkin(9) * t263;
t320 = -t211 * t278 + t255 * t237 + t323;
t212 = pkin(4) * t247 + pkin(9) * t246;
t319 = t278 * t212 - t237 * t256 + t322;
t318 = t256 * t211 - t255 * t212 + t321;
t306 = Icges(2,4) * t312;
t297 = rSges(2,1) * t312 - rSges(2,2) * t309;
t296 = rSges(2,1) * t309 + rSges(2,2) * t312;
t295 = Icges(2,1) * t312 - t352;
t294 = Icges(2,1) * t309 + t306;
t293 = -Icges(2,2) * t309 + t306;
t292 = Icges(2,2) * t312 + t352;
t287 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t286 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t285 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t265 = t313 * rSges(3,3) + (rSges(3,1) * t316 + rSges(3,2) * t317) * t310;
t262 = Icges(3,5) * t313 + (Icges(3,1) * t316 + Icges(3,4) * t317) * t310;
t260 = Icges(3,3) * t313 + (Icges(3,5) * t316 + Icges(3,6) * t317) * t310;
t259 = V_base(5) * rSges(2,3) - t296 * V_base(6) + t338;
t258 = t297 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t250 = t296 * V_base(4) - t297 * V_base(5) + t334;
t241 = qJD(5) * t263 + t278;
t240 = rSges(3,1) * t277 - rSges(3,2) * t276 + rSges(3,3) * t350;
t239 = rSges(3,1) * t275 - rSges(3,2) * t274 - rSges(3,3) * t349;
t236 = t273 * rSges(4,1) + t272 * rSges(4,2) - rSges(4,3) * t347;
t235 = Icges(3,1) * t277 - Icges(3,4) * t276 + Icges(3,5) * t350;
t234 = Icges(3,1) * t275 - Icges(3,4) * t274 - Icges(3,5) * t349;
t231 = Icges(3,5) * t277 - Icges(3,6) * t276 + Icges(3,3) * t350;
t230 = Icges(3,5) * t275 - Icges(3,6) * t274 - Icges(3,3) * t349;
t229 = Icges(4,1) * t273 + Icges(4,4) * t272 - Icges(4,5) * t347;
t228 = Icges(4,4) * t273 + Icges(4,2) * t272 - Icges(4,6) * t347;
t224 = t264 * rSges(5,1) - t263 * rSges(5,2) - rSges(5,3) * t347;
t223 = Icges(5,1) * t264 - Icges(5,4) * t263 - Icges(5,5) * t347;
t222 = Icges(5,4) * t264 - Icges(5,2) * t263 - Icges(5,6) * t347;
t221 = Icges(5,5) * t264 - Icges(5,6) * t263 - Icges(5,3) * t347;
t215 = qJD(5) * t246 + t256;
t214 = qJD(5) * t244 + t255;
t209 = rSges(4,1) * t254 + rSges(4,2) * t253 + rSges(4,3) * t276;
t208 = rSges(4,1) * t252 + rSges(4,2) * t251 + rSges(4,3) * t274;
t207 = Icges(4,1) * t254 + Icges(4,4) * t253 + Icges(4,5) * t276;
t206 = Icges(4,1) * t252 + Icges(4,4) * t251 + Icges(4,5) * t274;
t205 = Icges(4,4) * t254 + Icges(4,2) * t253 + Icges(4,6) * t276;
t204 = Icges(4,4) * t252 + Icges(4,2) * t251 + Icges(4,6) * t274;
t198 = rSges(5,1) * t247 - rSges(5,2) * t246 + rSges(5,3) * t276;
t197 = rSges(5,1) * t245 - rSges(5,2) * t244 + rSges(5,3) * t274;
t196 = Icges(5,1) * t247 - Icges(5,4) * t246 + Icges(5,5) * t276;
t195 = Icges(5,1) * t245 - Icges(5,4) * t244 + Icges(5,5) * t274;
t194 = Icges(5,4) * t247 - Icges(5,2) * t246 + Icges(5,6) * t276;
t193 = Icges(5,4) * t245 - Icges(5,2) * t244 + Icges(5,6) * t274;
t192 = Icges(5,5) * t247 - Icges(5,6) * t246 + Icges(5,3) * t276;
t191 = Icges(5,5) * t245 - Icges(5,6) * t244 + Icges(5,3) * t274;
t190 = rSges(6,1) * t249 - rSges(6,2) * t248 + rSges(6,3) * t263;
t179 = -t239 * t300 + t265 * t288 + t329;
t178 = t240 * t300 - t265 * t289 + t327;
t175 = t239 * t289 - t240 * t288 + t328;
t174 = rSges(6,1) * t219 - rSges(6,2) * t218 + rSges(6,3) * t246;
t172 = rSges(6,1) * t217 - rSges(6,2) * t216 + rSges(6,3) * t244;
t158 = t236 * t288 + (-t208 - t242) * t300 + t326;
t157 = t209 * t300 + (-t236 - t281) * t289 + t325;
t156 = t289 * t208 + (-t209 - t243) * t288 + t324;
t155 = -t197 * t278 + t224 * t255 + t323;
t154 = t198 * t278 - t224 * t256 + t322;
t153 = t256 * t197 - t255 * t198 + t321;
t152 = -t172 * t241 + t190 * t214 + t320;
t151 = t174 * t241 - t190 * t215 + t319;
t150 = t215 * t172 - t214 * t174 + t318;
t149 = qJD(6) * t218 + t214 * t341 - t241 * t343 + t320;
t148 = qJD(6) * t216 - t215 * t341 + t241 * t342 + t319;
t147 = qJD(6) * t248 - t214 * t342 + t215 * t343 + t318;
t1 = m(7) * (t147 ^ 2 + t148 ^ 2 + t149 ^ 2) / 0.2e1 + m(6) * (t150 ^ 2 + t151 ^ 2 + t152 ^ 2) / 0.2e1 + m(5) * (t153 ^ 2 + t154 ^ 2 + t155 ^ 2) / 0.2e1 + m(4) * (t156 ^ 2 + t157 ^ 2 + t158 ^ 2) / 0.2e1 + m(3) * (t175 ^ 2 + t178 ^ 2 + t179 ^ 2) / 0.2e1 + t278 * ((-t192 * t347 - t263 * t194 + t264 * t196) * t256 + (-t191 * t347 - t263 * t193 + t264 * t195) * t255 + (-t221 * t347 - t263 * t222 + t264 * t223) * t278) / 0.2e1 + t256 * ((t192 * t276 - t194 * t246 + t196 * t247) * t256 + (t191 * t276 - t193 * t246 + t195 * t247) * t255 + (t221 * t276 - t222 * t246 + t223 * t247) * t278) / 0.2e1 + t255 * ((t192 * t274 - t194 * t244 + t196 * t245) * t256 + (t191 * t274 - t193 * t244 + t195 * t245) * t255 + (t221 * t274 - t222 * t244 + t223 * t245) * t278) / 0.2e1 + m(1) * (t285 ^ 2 + t286 ^ 2 + t287 ^ 2) / 0.2e1 + m(2) * (t250 ^ 2 + t258 ^ 2 + t259 ^ 2) / 0.2e1 + ((t216 * t368 + t217 * t366 + t244 * t367) * t241 + (t216 * t373 + t217 * t369 + t244 * t371) * t215 + (t374 * t216 + t370 * t217 + t372 * t244) * t214) * t214 / 0.2e1 + ((t218 * t368 + t219 * t366 + t246 * t367) * t241 + (t373 * t218 + t369 * t219 + t371 * t246) * t215 + (t218 * t374 + t219 * t370 + t246 * t372) * t214) * t215 / 0.2e1 + ((t248 * t368 + t249 * t366 + t263 * t367) * t241 + (t248 * t373 + t249 * t369 + t263 * t371) * t215 + (t248 * t374 + t249 * t370 + t263 * t372) * t214) * t241 / 0.2e1 + ((t228 * t251 + t229 * t252 - t260 * t349 + t262 * t275 + t274 * t363) * t300 + (t205 * t251 + t207 * t252 - t231 * t349 + t235 * t275 + t274 * t364) * t289 + (t204 * t251 + t206 * t252 - t230 * t349 + t234 * t275 + t274 * t365) * t288) * t288 / 0.2e1 + ((t228 * t253 + t229 * t254 + t260 * t350 + t262 * t277 + t276 * t363) * t300 + (t205 * t253 + t207 * t254 + t231 * t350 + t235 * t277 + t276 * t364) * t289 + (t204 * t253 + t206 * t254 + t230 * t350 + t234 * t277 + t276 * t365) * t288) * t289 / 0.2e1 + ((t230 * t288 + t231 * t289 + t260 * t300) * t313 + ((t233 * t317 + t235 * t316) * t289 + (t232 * t317 + t234 * t316) * t288 + (t261 * t317 + t262 * t316) * t300) * t310 + (-t203 * t347 + t272 * t205 + t273 * t207) * t289 + (-t202 * t347 + t272 * t204 + t273 * t206) * t288 + (-t227 * t347 + t272 * t228 + t273 * t229) * t300) * t300 / 0.2e1 + ((-t292 * t309 + t294 * t312 + Icges(1,4)) * V_base(5) + (-t293 * t309 + t295 * t312 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t292 * t312 + t294 * t309 + Icges(1,2)) * V_base(5) + (t293 * t312 + t295 * t309 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t312 - Icges(2,6) * t309 + Icges(1,5)) * V_base(4) + (Icges(2,5) * t309 + Icges(2,6) * t312 + Icges(1,6)) * V_base(5) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
