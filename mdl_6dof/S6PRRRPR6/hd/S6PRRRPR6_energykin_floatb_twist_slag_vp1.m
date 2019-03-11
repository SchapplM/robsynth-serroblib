% Calculate kinetic energy for
% S6PRRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-03-08 23:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRPR6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR6_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRRPR6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR6_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPR6_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRPR6_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:31:03
% EndTime: 2019-03-08 23:31:09
% DurationCPUTime: 6.22s
% Computational Cost: add. (3189->392), mult. (7634->565), div. (0->0), fcn. (9432->12), ass. (0->170)
t365 = Icges(5,1) + Icges(6,1);
t364 = -Icges(5,4) + Icges(6,5);
t363 = Icges(6,4) + Icges(5,5);
t362 = Icges(5,2) + Icges(6,3);
t361 = Icges(6,2) + Icges(5,3);
t360 = -Icges(5,6) + Icges(6,6);
t307 = sin(pkin(11));
t309 = cos(pkin(11));
t315 = cos(qJ(2));
t313 = sin(qJ(2));
t343 = cos(pkin(6));
t329 = t313 * t343;
t274 = t307 * t315 + t309 * t329;
t308 = sin(pkin(6));
t312 = sin(qJ(3));
t339 = t308 * t312;
t345 = cos(qJ(3));
t256 = t274 * t345 - t309 * t339;
t328 = t315 * t343;
t273 = t307 * t313 - t309 * t328;
t311 = sin(qJ(4));
t344 = cos(qJ(4));
t225 = t256 * t311 - t273 * t344;
t226 = t256 * t344 + t273 * t311;
t331 = t308 * t345;
t255 = t274 * t312 + t309 * t331;
t357 = t362 * t225 + t364 * t226 + t360 * t255;
t276 = -t307 * t329 + t309 * t315;
t258 = t276 * t345 + t307 * t339;
t275 = t307 * t328 + t309 * t313;
t227 = t258 * t311 - t275 * t344;
t228 = t258 * t344 + t275 * t311;
t257 = t276 * t312 - t307 * t331;
t356 = t362 * t227 + t364 * t228 + t360 * t257;
t355 = t360 * t225 + t363 * t226 + t361 * t255;
t354 = t360 * t227 + t363 * t228 + t361 * t257;
t353 = t364 * t225 + t365 * t226 + t363 * t255;
t352 = t364 * t227 + t365 * t228 + t363 * t257;
t281 = t312 * t343 + t313 * t331;
t338 = t308 * t315;
t259 = t281 * t311 + t338 * t344;
t260 = t281 * t344 - t311 * t338;
t280 = t313 * t339 - t343 * t345;
t351 = t362 * t259 + t364 * t260 + t360 * t280;
t350 = t360 * t259 + t363 * t260 + t361 * t280;
t349 = t364 * t259 + t365 * t260 + t363 * t280;
t342 = Icges(2,4) * t307;
t341 = t307 * t308;
t340 = t308 * t309;
t337 = qJD(2) * t308;
t336 = V_base(5) * qJ(1) + V_base(1);
t332 = qJD(1) + V_base(3);
t330 = t343 * pkin(7);
t289 = t307 * t337 + V_base(4);
t300 = qJD(2) * t343 + V_base(6);
t254 = qJD(3) * t275 + t289;
t222 = qJD(4) * t257 + t254;
t288 = -t309 * t337 + V_base(5);
t253 = qJD(3) * t273 + t288;
t277 = -qJD(3) * t338 + t300;
t283 = pkin(1) * t307 - pkin(7) * t340;
t327 = -t283 * V_base(6) + V_base(5) * t330 + t336;
t284 = pkin(1) * t309 + pkin(7) * t341;
t326 = V_base(4) * t283 - t284 * V_base(5) + t332;
t221 = qJD(4) * t255 + t253;
t248 = qJD(4) * t280 + t277;
t325 = V_base(6) * t284 + V_base(2) + (-t330 - qJ(1)) * V_base(4);
t245 = pkin(2) * t274 + pkin(8) * t273;
t282 = (pkin(2) * t313 - pkin(8) * t315) * t308;
t324 = -t245 * t300 + t288 * t282 + t327;
t246 = pkin(2) * t276 + pkin(8) * t275;
t323 = t289 * t245 - t246 * t288 + t326;
t322 = t300 * t246 - t289 * t282 + t325;
t218 = pkin(3) * t256 + pkin(9) * t255;
t247 = pkin(3) * t281 + pkin(9) * t280;
t321 = -t218 * t277 + t253 * t247 + t324;
t219 = pkin(3) * t258 + pkin(9) * t257;
t320 = t254 * t218 - t219 * t253 + t323;
t220 = pkin(4) * t260 + qJ(5) * t259;
t319 = qJD(5) * t227 + t221 * t220 + t321;
t191 = pkin(4) * t226 + qJ(5) * t225;
t318 = qJD(5) * t259 + t222 * t191 + t320;
t317 = t277 * t219 - t254 * t247 + t322;
t192 = pkin(4) * t228 + qJ(5) * t227;
t316 = qJD(5) * t225 + t248 * t192 + t317;
t314 = cos(qJ(6));
t310 = sin(qJ(6));
t305 = Icges(2,4) * t309;
t297 = rSges(2,1) * t309 - rSges(2,2) * t307;
t296 = rSges(2,1) * t307 + rSges(2,2) * t309;
t295 = Icges(2,1) * t309 - t342;
t294 = Icges(2,1) * t307 + t305;
t293 = -Icges(2,2) * t307 + t305;
t292 = Icges(2,2) * t309 + t342;
t287 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t286 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t285 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t267 = t343 * rSges(3,3) + (rSges(3,1) * t313 + rSges(3,2) * t315) * t308;
t266 = Icges(3,5) * t343 + (Icges(3,1) * t313 + Icges(3,4) * t315) * t308;
t265 = Icges(3,6) * t343 + (Icges(3,4) * t313 + Icges(3,2) * t315) * t308;
t264 = Icges(3,3) * t343 + (Icges(3,5) * t313 + Icges(3,6) * t315) * t308;
t263 = V_base(5) * rSges(2,3) - t296 * V_base(6) + t336;
t262 = t297 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t251 = t296 * V_base(4) - t297 * V_base(5) + t332;
t244 = t281 * rSges(4,1) - t280 * rSges(4,2) - rSges(4,3) * t338;
t243 = Icges(4,1) * t281 - Icges(4,4) * t280 - Icges(4,5) * t338;
t242 = Icges(4,4) * t281 - Icges(4,2) * t280 - Icges(4,6) * t338;
t241 = Icges(4,5) * t281 - Icges(4,6) * t280 - Icges(4,3) * t338;
t240 = rSges(3,1) * t276 - rSges(3,2) * t275 + rSges(3,3) * t341;
t239 = rSges(3,1) * t274 - rSges(3,2) * t273 - rSges(3,3) * t340;
t238 = Icges(3,1) * t276 - Icges(3,4) * t275 + Icges(3,5) * t341;
t237 = Icges(3,1) * t274 - Icges(3,4) * t273 - Icges(3,5) * t340;
t236 = Icges(3,4) * t276 - Icges(3,2) * t275 + Icges(3,6) * t341;
t235 = Icges(3,4) * t274 - Icges(3,2) * t273 - Icges(3,6) * t340;
t234 = Icges(3,5) * t276 - Icges(3,6) * t275 + Icges(3,3) * t341;
t233 = Icges(3,5) * t274 - Icges(3,6) * t273 - Icges(3,3) * t340;
t230 = pkin(5) * t260 - pkin(10) * t280;
t229 = -qJD(6) * t280 + t248;
t217 = t259 * t310 + t260 * t314;
t216 = t259 * t314 - t260 * t310;
t214 = rSges(5,1) * t260 - rSges(5,2) * t259 + rSges(5,3) * t280;
t213 = rSges(6,1) * t260 + rSges(6,2) * t280 + rSges(6,3) * t259;
t205 = rSges(4,1) * t258 - rSges(4,2) * t257 + rSges(4,3) * t275;
t204 = rSges(4,1) * t256 - rSges(4,2) * t255 + rSges(4,3) * t273;
t203 = Icges(4,1) * t258 - Icges(4,4) * t257 + Icges(4,5) * t275;
t202 = Icges(4,1) * t256 - Icges(4,4) * t255 + Icges(4,5) * t273;
t201 = Icges(4,4) * t258 - Icges(4,2) * t257 + Icges(4,6) * t275;
t200 = Icges(4,4) * t256 - Icges(4,2) * t255 + Icges(4,6) * t273;
t199 = Icges(4,5) * t258 - Icges(4,6) * t257 + Icges(4,3) * t275;
t198 = Icges(4,5) * t256 - Icges(4,6) * t255 + Icges(4,3) * t273;
t197 = pkin(5) * t228 - pkin(10) * t257;
t196 = pkin(5) * t226 - pkin(10) * t255;
t195 = -qJD(6) * t257 + t222;
t194 = -qJD(6) * t255 + t221;
t190 = t227 * t310 + t228 * t314;
t189 = t227 * t314 - t228 * t310;
t188 = t225 * t310 + t226 * t314;
t187 = t225 * t314 - t226 * t310;
t186 = -t239 * t300 + t267 * t288 + t327;
t185 = t300 * t240 - t289 * t267 + t325;
t183 = rSges(5,1) * t228 - rSges(5,2) * t227 + rSges(5,3) * t257;
t182 = rSges(6,1) * t228 + rSges(6,2) * t257 + rSges(6,3) * t227;
t181 = rSges(5,1) * t226 - rSges(5,2) * t225 + rSges(5,3) * t255;
t180 = rSges(6,1) * t226 + rSges(6,2) * t255 + rSges(6,3) * t225;
t166 = t239 * t289 - t240 * t288 + t326;
t165 = rSges(7,1) * t217 + rSges(7,2) * t216 - rSges(7,3) * t280;
t164 = Icges(7,1) * t217 + Icges(7,4) * t216 - Icges(7,5) * t280;
t163 = Icges(7,4) * t217 + Icges(7,2) * t216 - Icges(7,6) * t280;
t162 = Icges(7,5) * t217 + Icges(7,6) * t216 - Icges(7,3) * t280;
t160 = rSges(7,1) * t190 + rSges(7,2) * t189 - rSges(7,3) * t257;
t159 = rSges(7,1) * t188 + rSges(7,2) * t187 - rSges(7,3) * t255;
t158 = Icges(7,1) * t190 + Icges(7,4) * t189 - Icges(7,5) * t257;
t157 = Icges(7,1) * t188 + Icges(7,4) * t187 - Icges(7,5) * t255;
t156 = Icges(7,4) * t190 + Icges(7,2) * t189 - Icges(7,6) * t257;
t155 = Icges(7,4) * t188 + Icges(7,2) * t187 - Icges(7,6) * t255;
t154 = Icges(7,5) * t190 + Icges(7,6) * t189 - Icges(7,3) * t257;
t153 = Icges(7,5) * t188 + Icges(7,6) * t187 - Icges(7,3) * t255;
t152 = -t204 * t277 + t244 * t253 + t324;
t151 = t277 * t205 - t254 * t244 + t322;
t150 = t204 * t254 - t205 * t253 + t323;
t149 = -t181 * t248 + t214 * t221 + t321;
t148 = t248 * t183 - t222 * t214 + t317;
t147 = t181 * t222 - t183 * t221 + t320;
t146 = t213 * t221 + (-t180 - t191) * t248 + t319;
t145 = t248 * t182 + (-t213 - t220) * t222 + t316;
t144 = t180 * t222 + (-t182 - t192) * t221 + t318;
t143 = -t159 * t229 + t165 * t194 + t221 * t230 + (-t191 - t196) * t248 + t319;
t142 = t229 * t160 - t195 * t165 + t248 * t197 + (-t220 - t230) * t222 + t316;
t141 = t159 * t195 - t160 * t194 + t196 * t222 + (-t192 - t197) * t221 + t318;
t1 = t229 * ((-t154 * t280 + t156 * t216 + t158 * t217) * t195 + (-t153 * t280 + t155 * t216 + t157 * t217) * t194 + (-t162 * t280 + t163 * t216 + t164 * t217) * t229) / 0.2e1 + m(1) * (t285 ^ 2 + t286 ^ 2 + t287 ^ 2) / 0.2e1 + t254 * ((t199 * t275 - t201 * t257 + t203 * t258) * t254 + (t198 * t275 - t200 * t257 + t202 * t258) * t253 + (t241 * t275 - t242 * t257 + t243 * t258) * t277) / 0.2e1 + t253 * ((t199 * t273 - t201 * t255 + t203 * t256) * t254 + (t198 * t273 - t200 * t255 + t202 * t256) * t253 + (t241 * t273 - t242 * t255 + t243 * t256) * t277) / 0.2e1 + m(2) * (t251 ^ 2 + t262 ^ 2 + t263 ^ 2) / 0.2e1 + t195 * ((-t257 * t154 + t189 * t156 + t190 * t158) * t195 + (-t153 * t257 + t155 * t189 + t157 * t190) * t194 + (-t162 * t257 + t163 * t189 + t164 * t190) * t229) / 0.2e1 + t194 * ((-t154 * t255 + t156 * t187 + t158 * t188) * t195 + (-t255 * t153 + t187 * t155 + t188 * t157) * t194 + (-t162 * t255 + t163 * t187 + t164 * t188) * t229) / 0.2e1 + m(3) * (t166 ^ 2 + t185 ^ 2 + t186 ^ 2) / 0.2e1 + m(5) * (t147 ^ 2 + t148 ^ 2 + t149 ^ 2) / 0.2e1 + m(4) * (t150 ^ 2 + t151 ^ 2 + t152 ^ 2) / 0.2e1 + m(6) * (t144 ^ 2 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + m(7) * (t141 ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + t288 * ((-t234 * t340 - t236 * t273 + t238 * t274) * t289 + (-t233 * t340 - t235 * t273 + t237 * t274) * t288 + (-t264 * t340 - t265 * t273 + t266 * t274) * t300) / 0.2e1 + t289 * ((t234 * t341 - t236 * t275 + t238 * t276) * t289 + (t233 * t341 - t235 * t275 + t237 * t276) * t288 + (t264 * t341 - t265 * t275 + t266 * t276) * t300) / 0.2e1 + t277 * ((-t199 * t338 - t280 * t201 + t281 * t203) * t254 + (-t198 * t338 - t280 * t200 + t281 * t202) * t253 + (-t241 * t338 - t280 * t242 + t281 * t243) * t277) / 0.2e1 + t300 * (((t236 * t315 + t238 * t313) * t289 + (t235 * t315 + t237 * t313) * t288 + (t265 * t315 + t266 * t313) * t300) * t308 + (t233 * t288 + t234 * t289 + t264 * t300) * t343) / 0.2e1 + ((t225 * t351 + t226 * t349 + t255 * t350) * t248 + (t225 * t356 + t226 * t352 + t255 * t354) * t222 + (t225 * t357 + t226 * t353 + t255 * t355) * t221) * t221 / 0.2e1 + ((t227 * t351 + t228 * t349 + t257 * t350) * t248 + (t227 * t356 + t228 * t352 + t257 * t354) * t222 + (t227 * t357 + t228 * t353 + t257 * t355) * t221) * t222 / 0.2e1 + ((t259 * t351 + t260 * t349 + t280 * t350) * t248 + (t259 * t356 + t260 * t352 + t280 * t354) * t222 + (t259 * t357 + t260 * t353 + t280 * t355) * t221) * t248 / 0.2e1 + ((-t292 * t307 + t294 * t309 + Icges(1,4)) * V_base(5) + (-t293 * t307 + t295 * t309 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t292 * t309 + t294 * t307 + Icges(1,2)) * V_base(5) + (t293 * t309 + t295 * t307 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t307 + Icges(2,6) * t309 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t309 - Icges(2,6) * t307 + Icges(1,5)) * V_base(4) + (Icges(2,3) / 0.2e1 + Icges(1,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
