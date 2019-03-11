% Calculate kinetic energy for
% S6PRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 21:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPRR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR1_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRPRR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR1_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR1_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPRR1_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:50:56
% EndTime: 2019-03-08 21:51:01
% DurationCPUTime: 5.00s
% Computational Cost: add. (3632->438), mult. (5494->632), div. (0->0), fcn. (6388->14), ass. (0->189)
t375 = Icges(4,3) + Icges(5,3);
t316 = sin(pkin(11));
t318 = cos(pkin(11));
t326 = cos(qJ(2));
t319 = cos(pkin(6));
t323 = sin(qJ(2));
t354 = t319 * t323;
t274 = t316 * t326 + t318 * t354;
t315 = qJ(3) + pkin(12);
t309 = sin(t315);
t310 = cos(t315);
t317 = sin(pkin(6));
t360 = t317 * t318;
t244 = -t274 * t309 - t310 * t360;
t245 = t274 * t310 - t309 * t360;
t322 = sin(qJ(3));
t325 = cos(qJ(3));
t357 = t317 * t325;
t251 = -t274 * t322 - t318 * t357;
t359 = t317 * t322;
t343 = t318 * t359;
t252 = t274 * t325 - t343;
t353 = t319 * t326;
t273 = t316 * t323 - t318 * t353;
t374 = Icges(4,5) * t252 + Icges(5,5) * t245 + Icges(4,6) * t251 + Icges(5,6) * t244 + t375 * t273;
t276 = -t316 * t354 + t318 * t326;
t361 = t316 * t317;
t246 = -t276 * t309 + t310 * t361;
t247 = t276 * t310 + t309 * t361;
t253 = -t276 * t322 + t316 * t357;
t344 = t316 * t359;
t254 = t276 * t325 + t344;
t275 = t316 * t353 + t318 * t323;
t373 = Icges(4,5) * t254 + Icges(5,5) * t247 + Icges(4,6) * t253 + Icges(5,6) * t246 + t375 * t275;
t358 = t317 * t323;
t264 = -t309 * t358 + t310 * t319;
t265 = t309 * t319 + t310 * t358;
t280 = t319 * t325 - t322 * t358;
t355 = t319 * t322;
t281 = t323 * t357 + t355;
t356 = t317 * t326;
t372 = Icges(4,5) * t281 + Icges(5,5) * t265 + Icges(4,6) * t280 + Icges(5,6) * t264 - t375 * t356;
t365 = pkin(7) * t319;
t364 = t325 * pkin(3);
t362 = Icges(2,4) * t316;
t352 = pkin(4) * t310;
t350 = qJD(2) * t317;
t349 = V_base(5) * qJ(1) + V_base(1);
t345 = qJD(1) + V_base(3);
t292 = t316 * t350 + V_base(4);
t303 = qJD(2) * t319 + V_base(6);
t342 = qJ(5) + t315;
t341 = pkin(4) * t309;
t250 = qJD(3) * t275 + t292;
t340 = cos(t342);
t221 = qJD(5) * t275 + t250;
t291 = -t318 * t350 + V_base(5);
t339 = t317 * t340;
t249 = qJD(3) * t273 + t291;
t283 = pkin(1) * t316 - pkin(7) * t360;
t338 = -t283 * V_base(6) + V_base(5) * t365 + t349;
t284 = pkin(1) * t318 + pkin(7) * t361;
t337 = V_base(4) * t283 - V_base(5) * t284 + t345;
t220 = qJD(5) * t273 + t249;
t258 = (-qJD(3) - qJD(5)) * t356 + t303;
t336 = V_base(6) * t284 + V_base(2) + (-qJ(1) - t365) * V_base(4);
t240 = pkin(2) * t274 + pkin(8) * t273;
t282 = (pkin(2) * t323 - pkin(8) * t326) * t317;
t335 = -t240 * t303 + t291 * t282 + t338;
t241 = pkin(2) * t276 + pkin(8) * t275;
t334 = t292 * t240 - t291 * t241 + t337;
t333 = t303 * t241 - t282 * t292 + t336;
t234 = pkin(3) * t355 + (-qJ(4) * t326 + t323 * t364) * t317;
t332 = qJD(4) * t275 + t249 * t234 + t335;
t190 = pkin(3) * t344 + qJ(4) * t275 + t276 * t364;
t277 = -qJD(3) * t356 + t303;
t331 = qJD(4) * t273 + t277 * t190 + t333;
t189 = -pkin(3) * t343 + qJ(4) * t273 + t274 * t364;
t330 = -qJD(4) * t356 + t250 * t189 + t334;
t162 = pkin(9) * t273 + t274 * t352 - t341 * t360;
t204 = t341 * t319 + (-pkin(9) * t326 + t323 * t352) * t317;
t329 = t249 * t204 + (-t162 - t189) * t277 + t332;
t163 = pkin(9) * t275 + t276 * t352 + t341 * t361;
t328 = t277 * t163 + (-t204 - t234) * t250 + t331;
t327 = t250 * t162 + (-t163 - t190) * t249 + t330;
t324 = cos(qJ(6));
t321 = sin(qJ(6));
t311 = Icges(2,4) * t318;
t306 = sin(t342);
t300 = rSges(2,1) * t318 - rSges(2,2) * t316;
t299 = rSges(2,1) * t316 + rSges(2,2) * t318;
t298 = Icges(2,1) * t318 - t362;
t297 = Icges(2,1) * t316 + t311;
t296 = -Icges(2,2) * t316 + t311;
t295 = Icges(2,2) * t318 + t362;
t289 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t288 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t287 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t266 = t319 * rSges(3,3) + (rSges(3,1) * t323 + rSges(3,2) * t326) * t317;
t263 = Icges(3,5) * t319 + (Icges(3,1) * t323 + Icges(3,4) * t326) * t317;
t262 = Icges(3,6) * t319 + (Icges(3,4) * t323 + Icges(3,2) * t326) * t317;
t261 = Icges(3,3) * t319 + (Icges(3,5) * t323 + Icges(3,6) * t326) * t317;
t260 = t319 * t306 + t323 * t339;
t259 = t306 * t358 - t319 * t340;
t257 = V_base(5) * rSges(2,3) - t299 * V_base(6) + t349;
t256 = t300 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t248 = t299 * V_base(4) - t300 * V_base(5) + t345;
t243 = t260 * t324 - t321 * t356;
t242 = -t260 * t321 - t324 * t356;
t239 = t276 * t340 + t306 * t361;
t238 = t276 * t306 - t316 * t339;
t237 = t274 * t340 - t306 * t360;
t236 = t274 * t306 + t318 * t339;
t235 = t281 * rSges(4,1) + t280 * rSges(4,2) - rSges(4,3) * t356;
t233 = Icges(4,1) * t281 + Icges(4,4) * t280 - Icges(4,5) * t356;
t232 = Icges(4,4) * t281 + Icges(4,2) * t280 - Icges(4,6) * t356;
t230 = rSges(3,1) * t276 - rSges(3,2) * t275 + rSges(3,3) * t361;
t229 = rSges(3,1) * t274 - rSges(3,2) * t273 - rSges(3,3) * t360;
t228 = Icges(3,1) * t276 - Icges(3,4) * t275 + Icges(3,5) * t361;
t227 = Icges(3,1) * t274 - Icges(3,4) * t273 - Icges(3,5) * t360;
t226 = Icges(3,4) * t276 - Icges(3,2) * t275 + Icges(3,6) * t361;
t225 = Icges(3,4) * t274 - Icges(3,2) * t273 - Icges(3,6) * t360;
t224 = Icges(3,5) * t276 - Icges(3,6) * t275 + Icges(3,3) * t361;
t223 = Icges(3,5) * t274 - Icges(3,6) * t273 - Icges(3,3) * t360;
t218 = qJD(6) * t259 + t258;
t217 = pkin(5) * t260 + pkin(10) * t259;
t216 = t265 * rSges(5,1) + t264 * rSges(5,2) - rSges(5,3) * t356;
t215 = Icges(5,1) * t265 + Icges(5,4) * t264 - Icges(5,5) * t356;
t214 = Icges(5,4) * t265 + Icges(5,2) * t264 - Icges(5,6) * t356;
t212 = t260 * rSges(6,1) - t259 * rSges(6,2) - rSges(6,3) * t356;
t211 = Icges(6,1) * t260 - Icges(6,4) * t259 - Icges(6,5) * t356;
t210 = Icges(6,4) * t260 - Icges(6,2) * t259 - Icges(6,6) * t356;
t209 = Icges(6,5) * t260 - Icges(6,6) * t259 - Icges(6,3) * t356;
t208 = t239 * t324 + t275 * t321;
t207 = -t239 * t321 + t275 * t324;
t206 = t237 * t324 + t273 * t321;
t205 = -t237 * t321 + t273 * t324;
t202 = pkin(5) * t239 + pkin(10) * t238;
t201 = pkin(5) * t237 + pkin(10) * t236;
t200 = rSges(4,1) * t254 + rSges(4,2) * t253 + rSges(4,3) * t275;
t199 = rSges(4,1) * t252 + rSges(4,2) * t251 + rSges(4,3) * t273;
t198 = Icges(4,1) * t254 + Icges(4,4) * t253 + Icges(4,5) * t275;
t197 = Icges(4,1) * t252 + Icges(4,4) * t251 + Icges(4,5) * t273;
t196 = Icges(4,4) * t254 + Icges(4,2) * t253 + Icges(4,6) * t275;
t195 = Icges(4,4) * t252 + Icges(4,2) * t251 + Icges(4,6) * t273;
t192 = qJD(6) * t238 + t221;
t191 = qJD(6) * t236 + t220;
t188 = rSges(5,1) * t247 + rSges(5,2) * t246 + rSges(5,3) * t275;
t187 = rSges(5,1) * t245 + rSges(5,2) * t244 + rSges(5,3) * t273;
t186 = Icges(5,1) * t247 + Icges(5,4) * t246 + Icges(5,5) * t275;
t185 = Icges(5,1) * t245 + Icges(5,4) * t244 + Icges(5,5) * t273;
t184 = Icges(5,4) * t247 + Icges(5,2) * t246 + Icges(5,6) * t275;
t183 = Icges(5,4) * t245 + Icges(5,2) * t244 + Icges(5,6) * t273;
t180 = rSges(6,1) * t239 - rSges(6,2) * t238 + rSges(6,3) * t275;
t179 = rSges(6,1) * t237 - rSges(6,2) * t236 + rSges(6,3) * t273;
t178 = Icges(6,1) * t239 - Icges(6,4) * t238 + Icges(6,5) * t275;
t177 = Icges(6,1) * t237 - Icges(6,4) * t236 + Icges(6,5) * t273;
t176 = Icges(6,4) * t239 - Icges(6,2) * t238 + Icges(6,6) * t275;
t175 = Icges(6,4) * t237 - Icges(6,2) * t236 + Icges(6,6) * t273;
t174 = Icges(6,5) * t239 - Icges(6,6) * t238 + Icges(6,3) * t275;
t173 = Icges(6,5) * t237 - Icges(6,6) * t236 + Icges(6,3) * t273;
t172 = rSges(7,1) * t243 + rSges(7,2) * t242 + rSges(7,3) * t259;
t171 = Icges(7,1) * t243 + Icges(7,4) * t242 + Icges(7,5) * t259;
t170 = Icges(7,4) * t243 + Icges(7,2) * t242 + Icges(7,6) * t259;
t169 = Icges(7,5) * t243 + Icges(7,6) * t242 + Icges(7,3) * t259;
t166 = -t229 * t303 + t266 * t291 + t338;
t165 = t230 * t303 - t266 * t292 + t336;
t160 = t229 * t292 - t230 * t291 + t337;
t158 = rSges(7,1) * t208 + rSges(7,2) * t207 + rSges(7,3) * t238;
t157 = rSges(7,1) * t206 + rSges(7,2) * t205 + rSges(7,3) * t236;
t156 = Icges(7,1) * t208 + Icges(7,4) * t207 + Icges(7,5) * t238;
t155 = Icges(7,1) * t206 + Icges(7,4) * t205 + Icges(7,5) * t236;
t154 = Icges(7,4) * t208 + Icges(7,2) * t207 + Icges(7,6) * t238;
t153 = Icges(7,4) * t206 + Icges(7,2) * t205 + Icges(7,6) * t236;
t152 = Icges(7,5) * t208 + Icges(7,6) * t207 + Icges(7,3) * t238;
t151 = Icges(7,5) * t206 + Icges(7,6) * t205 + Icges(7,3) * t236;
t150 = -t199 * t277 + t235 * t249 + t335;
t149 = t200 * t277 - t235 * t250 + t333;
t148 = t199 * t250 - t200 * t249 + t334;
t147 = t216 * t249 + (-t187 - t189) * t277 + t332;
t146 = t188 * t277 + (-t216 - t234) * t250 + t331;
t145 = t250 * t187 + (-t188 - t190) * t249 + t330;
t144 = -t179 * t258 + t212 * t220 + t329;
t143 = t180 * t258 - t212 * t221 + t328;
t142 = t221 * t179 - t220 * t180 + t327;
t141 = -t157 * t218 + t172 * t191 - t201 * t258 + t217 * t220 + t329;
t140 = t158 * t218 - t172 * t192 + t202 * t258 - t217 * t221 + t328;
t139 = t192 * t157 - t191 * t158 + t221 * t201 - t220 * t202 + t327;
t1 = t303 * ((t223 * t291 + t224 * t292 + t261 * t303) * t319 + ((t226 * t326 + t228 * t323) * t292 + (t225 * t326 + t227 * t323) * t291 + (t262 * t326 + t263 * t323) * t303) * t317) / 0.2e1 + t291 * ((-t224 * t360 - t226 * t273 + t228 * t274) * t292 + (-t223 * t360 - t225 * t273 + t227 * t274) * t291 + (-t261 * t360 - t262 * t273 + t263 * t274) * t303) / 0.2e1 + t292 * ((t224 * t361 - t226 * t275 + t228 * t276) * t292 + (t223 * t361 - t225 * t275 + t227 * t276) * t291 + (t261 * t361 - t262 * t275 + t263 * t276) * t303) / 0.2e1 + t258 * ((-t174 * t356 - t259 * t176 + t260 * t178) * t221 + (-t173 * t356 - t259 * t175 + t260 * t177) * t220 + (-t209 * t356 - t259 * t210 + t260 * t211) * t258) / 0.2e1 + m(3) * (t160 ^ 2 + t165 ^ 2 + t166 ^ 2) / 0.2e1 + t191 * ((t152 * t236 + t154 * t205 + t156 * t206) * t192 + (t236 * t151 + t205 * t153 + t206 * t155) * t191 + (t169 * t236 + t170 * t205 + t171 * t206) * t218) / 0.2e1 + t192 * ((t238 * t152 + t207 * t154 + t208 * t156) * t192 + (t151 * t238 + t153 * t207 + t155 * t208) * t191 + (t169 * t238 + t170 * t207 + t171 * t208) * t218) / 0.2e1 + m(2) * (t248 ^ 2 + t256 ^ 2 + t257 ^ 2) / 0.2e1 + t218 * ((t152 * t259 + t154 * t242 + t156 * t243) * t192 + (t151 * t259 + t153 * t242 + t155 * t243) * t191 + (t259 * t169 + t242 * t170 + t243 * t171) * t218) / 0.2e1 + t220 * ((t174 * t273 - t176 * t236 + t178 * t237) * t221 + (t273 * t173 - t236 * t175 + t237 * t177) * t220 + (t209 * t273 - t210 * t236 + t211 * t237) * t258) / 0.2e1 + t221 * ((t275 * t174 - t238 * t176 + t239 * t178) * t221 + (t173 * t275 - t175 * t238 + t177 * t239) * t220 + (t209 * t275 - t210 * t238 + t211 * t239) * t258) / 0.2e1 + m(1) * (t287 ^ 2 + t288 ^ 2 + t289 ^ 2) / 0.2e1 + m(7) * (t139 ^ 2 + t140 ^ 2 + t141 ^ 2) / 0.2e1 + m(6) * (t142 ^ 2 + t143 ^ 2 + t144 ^ 2) / 0.2e1 + m(5) * (t145 ^ 2 + t146 ^ 2 + t147 ^ 2) / 0.2e1 + m(4) * (t148 ^ 2 + t149 ^ 2 + t150 ^ 2) / 0.2e1 + ((t214 * t244 + t215 * t245 + t232 * t251 + t233 * t252 + t273 * t372) * t277 + (t184 * t244 + t186 * t245 + t196 * t251 + t198 * t252 + t273 * t373) * t250 + (t183 * t244 + t185 * t245 + t195 * t251 + t197 * t252 + t374 * t273) * t249) * t249 / 0.2e1 + ((t214 * t246 + t215 * t247 + t232 * t253 + t233 * t254 + t275 * t372) * t277 + (t184 * t246 + t186 * t247 + t196 * t253 + t198 * t254 + t373 * t275) * t250 + (t183 * t246 + t185 * t247 + t195 * t253 + t197 * t254 + t275 * t374) * t249) * t250 / 0.2e1 + ((t264 * t214 + t265 * t215 + t280 * t232 + t281 * t233 - t356 * t372) * t277 + (t264 * t184 + t265 * t186 + t280 * t196 + t281 * t198 - t356 * t373) * t250 + (t264 * t183 + t265 * t185 + t280 * t195 + t281 * t197 - t356 * t374) * t249) * t277 / 0.2e1 + ((-t295 * t316 + t297 * t318 + Icges(1,4)) * V_base(5) + (-t296 * t316 + t298 * t318 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t295 * t318 + t297 * t316 + Icges(1,2)) * V_base(5) + (t296 * t318 + t298 * t316 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t318 - Icges(2,6) * t316 + Icges(1,5)) * V_base(4) + (Icges(2,5) * t316 + Icges(2,6) * t318 + Icges(1,6)) * V_base(5) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
