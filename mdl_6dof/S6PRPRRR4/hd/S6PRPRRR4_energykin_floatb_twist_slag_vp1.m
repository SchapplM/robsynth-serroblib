% Calculate kinetic energy for
% S6PRPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRRR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR4_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRPRRR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR4_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRR4_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRR4_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:35:50
% EndTime: 2019-03-08 20:35:54
% DurationCPUTime: 4.36s
% Computational Cost: add. (3635->442), mult. (5948->637), div. (0->0), fcn. (7076->14), ass. (0->194)
t307 = sin(pkin(11));
t310 = cos(pkin(11));
t316 = cos(qJ(2));
t311 = cos(pkin(6));
t314 = sin(qJ(2));
t344 = t311 * t314;
t269 = t307 * t316 + t310 * t344;
t306 = sin(pkin(12));
t309 = cos(pkin(12));
t308 = sin(pkin(6));
t347 = t308 * t310;
t246 = -t269 * t306 - t309 * t347;
t333 = t306 * t347;
t247 = t269 * t309 - t333;
t343 = t311 * t316;
t268 = t307 * t314 - t310 * t343;
t189 = Icges(4,5) * t247 + Icges(4,6) * t246 + Icges(4,3) * t268;
t223 = Icges(3,4) * t269 - Icges(3,2) * t268 - Icges(3,6) * t347;
t364 = t189 - t223;
t271 = -t307 * t344 + t310 * t316;
t348 = t307 * t308;
t248 = -t271 * t306 + t309 * t348;
t334 = t306 * t348;
t249 = t271 * t309 + t334;
t270 = t307 * t343 + t310 * t314;
t190 = Icges(4,5) * t249 + Icges(4,6) * t248 + Icges(4,3) * t270;
t224 = Icges(3,4) * t271 - Icges(3,2) * t270 + Icges(3,6) * t348;
t363 = t190 - t224;
t346 = t308 * t314;
t266 = -t306 * t346 + t309 * t311;
t349 = t306 * t311;
t267 = t309 * t346 + t349;
t345 = t308 * t316;
t218 = Icges(4,5) * t267 + Icges(4,6) * t266 - Icges(4,3) * t345;
t257 = Icges(3,6) * t311 + (Icges(3,4) * t314 + Icges(3,2) * t316) * t308;
t362 = t218 - t257;
t356 = pkin(7) * t311;
t355 = pkin(3) * t309;
t315 = cos(qJ(5));
t354 = pkin(5) * t315;
t352 = Icges(2,4) * t307;
t313 = sin(qJ(5));
t351 = t268 * t313;
t350 = t270 * t313;
t341 = qJD(2) * t308;
t340 = pkin(12) + qJ(4);
t339 = V_base(5) * qJ(1) + V_base(1);
t335 = qJD(1) + V_base(3);
t332 = t313 * t345;
t283 = t307 * t341 + V_base(4);
t294 = qJD(2) * t311 + V_base(6);
t331 = cos(t340);
t251 = qJD(4) * t270 + t283;
t330 = t308 * t331;
t300 = sin(t340);
t241 = t271 * t300 - t307 * t330;
t201 = qJD(5) * t241 + t251;
t282 = -t310 * t341 + V_base(5);
t250 = qJD(4) * t268 + t282;
t272 = -qJD(4) * t345 + t294;
t276 = pkin(1) * t307 - pkin(7) * t347;
t329 = -t276 * V_base(6) + V_base(5) * t356 + t339;
t277 = pkin(1) * t310 + pkin(7) * t348;
t328 = V_base(4) * t276 - t277 * V_base(5) + t335;
t239 = t269 * t300 + t310 * t330;
t200 = qJD(5) * t239 + t250;
t259 = t300 * t346 - t311 * t331;
t232 = qJD(5) * t259 + t272;
t327 = V_base(6) * t277 + V_base(2) + (-qJ(1) - t356) * V_base(4);
t275 = (pkin(2) * t314 - qJ(3) * t316) * t308;
t326 = qJD(3) * t270 + t282 * t275 + t329;
t236 = pkin(2) * t271 + qJ(3) * t270;
t325 = qJD(3) * t268 + t294 * t236 + t327;
t235 = pkin(2) * t269 + qJ(3) * t268;
t324 = -qJD(3) * t345 + t283 * t235 + t328;
t186 = -pkin(3) * t333 + pkin(8) * t268 + t269 * t355;
t229 = pkin(3) * t349 + (-pkin(8) * t316 + t314 * t355) * t308;
t323 = t282 * t229 + (-t186 - t235) * t294 + t326;
t187 = pkin(3) * t334 + pkin(8) * t270 + t271 * t355;
t322 = t294 * t187 + (-t229 - t275) * t283 + t325;
t321 = t283 * t186 + (-t187 - t236) * t282 + t324;
t240 = t269 * t331 - t300 * t347;
t198 = t240 * pkin(4) + t239 * pkin(9);
t260 = t311 * t300 + t314 * t330;
t228 = t260 * pkin(4) + t259 * pkin(9);
t320 = -t198 * t272 + t250 * t228 + t323;
t242 = t271 * t331 + t300 * t348;
t199 = t242 * pkin(4) + t241 * pkin(9);
t319 = t272 * t199 - t228 * t251 + t322;
t318 = t251 * t198 - t199 * t250 + t321;
t305 = qJ(5) + qJ(6);
t303 = cos(t305);
t302 = sin(t305);
t301 = Icges(2,4) * t310;
t291 = rSges(2,1) * t310 - rSges(2,2) * t307;
t290 = rSges(2,1) * t307 + rSges(2,2) * t310;
t289 = Icges(2,1) * t310 - t352;
t288 = Icges(2,1) * t307 + t301;
t287 = -Icges(2,2) * t307 + t301;
t286 = Icges(2,2) * t310 + t352;
t281 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t280 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t279 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t261 = rSges(3,3) * t311 + (rSges(3,1) * t314 + rSges(3,2) * t316) * t308;
t258 = Icges(3,5) * t311 + (Icges(3,1) * t314 + Icges(3,4) * t316) * t308;
t256 = Icges(3,3) * t311 + (Icges(3,5) * t314 + Icges(3,6) * t316) * t308;
t254 = V_base(5) * rSges(2,3) - t290 * V_base(6) + t339;
t253 = t291 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t245 = t290 * V_base(4) - t291 * V_base(5) + t335;
t244 = t260 * t315 - t332;
t243 = -t260 * t313 - t315 * t345;
t234 = t260 * t303 - t302 * t345;
t233 = -t260 * t302 - t303 * t345;
t231 = rSges(3,1) * t271 - rSges(3,2) * t270 + rSges(3,3) * t348;
t230 = rSges(3,1) * t269 - rSges(3,2) * t268 - rSges(3,3) * t347;
t227 = rSges(4,1) * t267 + rSges(4,2) * t266 - rSges(4,3) * t345;
t226 = Icges(3,1) * t271 - Icges(3,4) * t270 + Icges(3,5) * t348;
t225 = Icges(3,1) * t269 - Icges(3,4) * t268 - Icges(3,5) * t347;
t222 = Icges(3,5) * t271 - Icges(3,6) * t270 + Icges(3,3) * t348;
t221 = Icges(3,5) * t269 - Icges(3,6) * t268 - Icges(3,3) * t347;
t220 = Icges(4,1) * t267 + Icges(4,4) * t266 - Icges(4,5) * t345;
t219 = Icges(4,4) * t267 + Icges(4,2) * t266 - Icges(4,6) * t345;
t215 = rSges(5,1) * t260 - rSges(5,2) * t259 - rSges(5,3) * t345;
t214 = Icges(5,1) * t260 - Icges(5,4) * t259 - Icges(5,5) * t345;
t213 = Icges(5,4) * t260 - Icges(5,2) * t259 - Icges(5,6) * t345;
t212 = Icges(5,5) * t260 - Icges(5,6) * t259 - Icges(5,3) * t345;
t210 = t242 * t315 + t350;
t209 = -t242 * t313 + t270 * t315;
t208 = t240 * t315 + t351;
t207 = -t240 * t313 + t268 * t315;
t206 = qJD(6) * t259 + t232;
t205 = t242 * t303 + t270 * t302;
t204 = -t242 * t302 + t270 * t303;
t203 = t240 * t303 + t268 * t302;
t202 = -t240 * t302 + t268 * t303;
t196 = rSges(4,1) * t249 + rSges(4,2) * t248 + rSges(4,3) * t270;
t195 = rSges(4,1) * t247 + rSges(4,2) * t246 + rSges(4,3) * t268;
t194 = Icges(4,1) * t249 + Icges(4,4) * t248 + Icges(4,5) * t270;
t193 = Icges(4,1) * t247 + Icges(4,4) * t246 + Icges(4,5) * t268;
t192 = Icges(4,4) * t249 + Icges(4,2) * t248 + Icges(4,6) * t270;
t191 = Icges(4,4) * t247 + Icges(4,2) * t246 + Icges(4,6) * t268;
t185 = rSges(5,1) * t242 - rSges(5,2) * t241 + rSges(5,3) * t270;
t184 = rSges(5,1) * t240 - rSges(5,2) * t239 + rSges(5,3) * t268;
t183 = Icges(5,1) * t242 - Icges(5,4) * t241 + Icges(5,5) * t270;
t182 = Icges(5,1) * t240 - Icges(5,4) * t239 + Icges(5,5) * t268;
t181 = Icges(5,4) * t242 - Icges(5,2) * t241 + Icges(5,6) * t270;
t180 = Icges(5,4) * t240 - Icges(5,2) * t239 + Icges(5,6) * t268;
t179 = Icges(5,5) * t242 - Icges(5,6) * t241 + Icges(5,3) * t270;
t178 = Icges(5,5) * t240 - Icges(5,6) * t239 + Icges(5,3) * t268;
t177 = rSges(6,1) * t244 + rSges(6,2) * t243 + rSges(6,3) * t259;
t176 = Icges(6,1) * t244 + Icges(6,4) * t243 + Icges(6,5) * t259;
t175 = Icges(6,4) * t244 + Icges(6,2) * t243 + Icges(6,6) * t259;
t174 = Icges(6,5) * t244 + Icges(6,6) * t243 + Icges(6,3) * t259;
t172 = rSges(7,1) * t234 + rSges(7,2) * t233 + rSges(7,3) * t259;
t171 = Icges(7,1) * t234 + Icges(7,4) * t233 + Icges(7,5) * t259;
t170 = Icges(7,4) * t234 + Icges(7,2) * t233 + Icges(7,6) * t259;
t169 = Icges(7,5) * t234 + Icges(7,6) * t233 + Icges(7,3) * t259;
t167 = qJD(6) * t241 + t201;
t166 = qJD(6) * t239 + t200;
t164 = -pkin(5) * t332 + pkin(10) * t259 + t260 * t354;
t163 = -t230 * t294 + t261 * t282 + t329;
t162 = t231 * t294 - t261 * t283 + t327;
t161 = t230 * t283 - t231 * t282 + t328;
t160 = rSges(6,1) * t210 + rSges(6,2) * t209 + rSges(6,3) * t241;
t159 = rSges(6,1) * t208 + rSges(6,2) * t207 + rSges(6,3) * t239;
t158 = Icges(6,1) * t210 + Icges(6,4) * t209 + Icges(6,5) * t241;
t157 = Icges(6,1) * t208 + Icges(6,4) * t207 + Icges(6,5) * t239;
t156 = Icges(6,4) * t210 + Icges(6,2) * t209 + Icges(6,6) * t241;
t155 = Icges(6,4) * t208 + Icges(6,2) * t207 + Icges(6,6) * t239;
t154 = Icges(6,5) * t210 + Icges(6,6) * t209 + Icges(6,3) * t241;
t153 = Icges(6,5) * t208 + Icges(6,6) * t207 + Icges(6,3) * t239;
t152 = rSges(7,1) * t205 + rSges(7,2) * t204 + rSges(7,3) * t241;
t151 = rSges(7,1) * t203 + rSges(7,2) * t202 + rSges(7,3) * t239;
t150 = Icges(7,1) * t205 + Icges(7,4) * t204 + Icges(7,5) * t241;
t149 = Icges(7,1) * t203 + Icges(7,4) * t202 + Icges(7,5) * t239;
t148 = Icges(7,4) * t205 + Icges(7,2) * t204 + Icges(7,6) * t241;
t147 = Icges(7,4) * t203 + Icges(7,2) * t202 + Icges(7,6) * t239;
t146 = Icges(7,5) * t205 + Icges(7,6) * t204 + Icges(7,3) * t241;
t145 = Icges(7,5) * t203 + Icges(7,6) * t202 + Icges(7,3) * t239;
t144 = pkin(5) * t350 + pkin(10) * t241 + t242 * t354;
t143 = pkin(5) * t351 + pkin(10) * t239 + t240 * t354;
t142 = t227 * t282 + (-t195 - t235) * t294 + t326;
t141 = t196 * t294 + (-t227 - t275) * t283 + t325;
t140 = t195 * t283 + (-t196 - t236) * t282 + t324;
t139 = -t184 * t272 + t215 * t250 + t323;
t138 = t185 * t272 - t215 * t251 + t322;
t137 = t184 * t251 - t185 * t250 + t321;
t136 = -t159 * t232 + t177 * t200 + t320;
t135 = t160 * t232 - t177 * t201 + t319;
t134 = t159 * t201 - t160 * t200 + t318;
t133 = -t143 * t232 - t151 * t206 + t164 * t200 + t166 * t172 + t320;
t132 = t144 * t232 + t152 * t206 - t164 * t201 - t167 * t172 + t319;
t131 = t143 * t201 - t144 * t200 + t151 * t167 - t152 * t166 + t318;
t1 = m(3) * (t161 ^ 2 + t162 ^ 2 + t163 ^ 2) / 0.2e1 + m(5) * (t137 ^ 2 + t138 ^ 2 + t139 ^ 2) / 0.2e1 + m(4) * (t140 ^ 2 + t141 ^ 2 + t142 ^ 2) / 0.2e1 + m(6) * (t134 ^ 2 + t135 ^ 2 + t136 ^ 2) / 0.2e1 + m(7) * (t131 ^ 2 + t132 ^ 2 + t133 ^ 2) / 0.2e1 + m(1) * (t279 ^ 2 + t280 ^ 2 + t281 ^ 2) / 0.2e1 + t251 * ((t179 * t270 - t181 * t241 + t183 * t242) * t251 + (t178 * t270 - t180 * t241 + t182 * t242) * t250 + (t212 * t270 - t213 * t241 + t214 * t242) * t272) / 0.2e1 + t250 * ((t179 * t268 - t181 * t239 + t183 * t240) * t251 + (t178 * t268 - t180 * t239 + t182 * t240) * t250 + (t212 * t268 - t213 * t239 + t214 * t240) * t272) / 0.2e1 + t232 * ((t154 * t259 + t156 * t243 + t158 * t244) * t201 + (t153 * t259 + t155 * t243 + t157 * t244) * t200 + (t259 * t174 + t175 * t243 + t244 * t176) * t232) / 0.2e1 + t206 * ((t146 * t259 + t148 * t233 + t150 * t234) * t167 + (t145 * t259 + t147 * t233 + t149 * t234) * t166 + (t259 * t169 + t233 * t170 + t234 * t171) * t206) / 0.2e1 + m(2) * (t245 ^ 2 + t253 ^ 2 + t254 ^ 2) / 0.2e1 + t167 * ((t241 * t146 + t204 * t148 + t150 * t205) * t167 + (t145 * t241 + t147 * t204 + t149 * t205) * t166 + (t169 * t241 + t170 * t204 + t171 * t205) * t206) / 0.2e1 + t201 * ((t241 * t154 + t209 * t156 + t158 * t210) * t201 + (t153 * t241 + t155 * t209 + t157 * t210) * t200 + (t174 * t241 + t175 * t209 + t176 * t210) * t232) / 0.2e1 + t166 * ((t146 * t239 + t148 * t202 + t150 * t203) * t167 + (t145 * t239 + t202 * t147 + t203 * t149) * t166 + (t169 * t239 + t170 * t202 + t171 * t203) * t206) / 0.2e1 + t200 * ((t154 * t239 + t156 * t207 + t158 * t208) * t201 + (t239 * t153 + t155 * t207 + t157 * t208) * t200 + (t174 * t239 + t175 * t207 + t176 * t208) * t232) / 0.2e1 + t272 * ((-t179 * t345 - t181 * t259 + t183 * t260) * t251 + (-t178 * t345 - t180 * t259 + t182 * t260) * t250 + (-t212 * t345 - t213 * t259 + t214 * t260) * t272) / 0.2e1 + ((t219 * t246 + t220 * t247 - t256 * t347 + t258 * t269 + t268 * t362) * t294 + (t192 * t246 + t194 * t247 - t222 * t347 + t226 * t269 + t268 * t363) * t283 + (t191 * t246 + t193 * t247 - t221 * t347 + t225 * t269 + t364 * t268) * t282) * t282 / 0.2e1 + ((t219 * t248 + t220 * t249 + t256 * t348 + t258 * t271 + t270 * t362) * t294 + (t192 * t248 + t194 * t249 + t222 * t348 + t226 * t271 + t363 * t270) * t283 + (t191 * t248 + t193 * t249 + t221 * t348 + t225 * t271 + t270 * t364) * t282) * t283 / 0.2e1 + ((t221 * t282 + t222 * t283 + t256 * t294) * t311 + ((t224 * t316 + t226 * t314) * t283 + (t223 * t316 + t225 * t314) * t282 + (t257 * t316 + t258 * t314) * t294) * t308 + (-t190 * t345 + t192 * t266 + t194 * t267) * t283 + (-t189 * t345 + t191 * t266 + t193 * t267) * t282 + (-t218 * t345 + t219 * t266 + t220 * t267) * t294) * t294 / 0.2e1 + ((-t286 * t307 + t288 * t310 + Icges(1,4)) * V_base(5) + (-t287 * t307 + t289 * t310 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t286 * t310 + t288 * t307 + Icges(1,2)) * V_base(5) + (t287 * t310 + t289 * t307 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t310 - Icges(2,6) * t307 + Icges(1,5)) * V_base(4) + (Icges(2,5) * t307 + Icges(2,6) * t310 + Icges(1,6)) * V_base(5) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
