% Calculate kinetic energy for
% S6RRPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
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
% Datum: 2019-03-09 08:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPPR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR4_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPPPR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR4_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPPR4_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPPR4_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:17:11
% EndTime: 2019-03-09 08:17:14
% DurationCPUTime: 4.03s
% Computational Cost: add. (1236->323), mult. (2497->435), div. (0->0), fcn. (2555->8), ass. (0->156)
t334 = Icges(3,4) + Icges(4,6);
t333 = Icges(3,1) + Icges(4,2);
t332 = Icges(3,2) + Icges(4,3);
t248 = cos(qJ(2));
t331 = t334 * t248;
t245 = sin(qJ(2));
t330 = t334 * t245;
t329 = Icges(4,4) - Icges(3,5);
t328 = Icges(4,5) - Icges(3,6);
t327 = t332 * t245 - t331;
t326 = t333 * t248 - t330;
t325 = Icges(4,1) + Icges(3,3);
t324 = Icges(5,1) + Icges(6,1);
t323 = Icges(5,4) - Icges(6,5);
t322 = Icges(6,4) + Icges(5,5);
t321 = Icges(5,2) + Icges(6,3);
t320 = Icges(6,6) - Icges(5,6);
t319 = Icges(5,3) + Icges(6,2);
t246 = sin(qJ(1));
t249 = cos(qJ(1));
t318 = t327 * t246 - t328 * t249;
t317 = t328 * t246 + t327 * t249;
t316 = -t329 * t246 + t326 * t249;
t315 = t326 * t246 + t329 * t249;
t314 = -t332 * t248 - t330;
t313 = t333 * t245 + t331;
t312 = t328 * t245 - t329 * t248;
t243 = cos(pkin(9));
t290 = t245 * t249;
t242 = sin(pkin(9));
t291 = t242 * t246;
t191 = -t243 * t290 + t291;
t289 = t246 * t243;
t192 = t242 * t290 + t289;
t287 = t248 * t249;
t311 = t191 * t321 - t192 * t323 + t287 * t320;
t193 = t242 * t249 + t245 * t289;
t194 = -t243 * t249 + t245 * t291;
t288 = t246 * t248;
t310 = -t193 * t321 - t194 * t323 + t288 * t320;
t309 = t191 * t320 + t192 * t322 + t287 * t319;
t308 = -t193 * t320 + t194 * t322 + t288 * t319;
t307 = -t191 * t323 + t192 * t324 + t287 * t322;
t306 = t193 * t323 + t194 * t324 + t288 * t322;
t305 = (t242 * t323 + t243 * t321) * t248 + t320 * t245;
t304 = (-t242 * t322 + t243 * t320) * t248 + t319 * t245;
t303 = (-t242 * t324 - t243 * t323) * t248 + t322 * t245;
t230 = -qJD(2) * t249 + V_base(5);
t231 = qJD(2) * t246 + V_base(4);
t237 = V_base(6) + qJD(1);
t302 = (t245 * t314 + t248 * t313) * t237 + (t245 * t317 + t248 * t316) * t231 + (t245 * t318 + t248 * t315) * t230;
t301 = (-t329 * t245 - t328 * t248) * t237 + (t246 * t325 + t312 * t249) * t231 + (t312 * t246 - t249 * t325) * t230;
t297 = Icges(2,4) * t246;
t292 = qJ(4) * t245;
t270 = pkin(2) * t248 + qJ(3) * t245;
t198 = t270 * t246;
t227 = t246 * pkin(1) - pkin(7) * t249;
t286 = -t198 - t227;
t199 = t270 * t249;
t204 = t246 * pkin(3) + qJ(4) * t287;
t285 = -t199 - t204;
t284 = qJD(3) * t245;
t283 = qJD(4) * t248;
t282 = qJD(6) * t248;
t281 = V_base(5) * pkin(6) + V_base(1);
t150 = pkin(4) * t192 + qJ(5) * t191;
t278 = -t150 + t285;
t205 = -pkin(3) * t249 + qJ(4) * t288;
t277 = -t205 + t286;
t222 = pkin(2) * t245 - qJ(3) * t248;
t276 = -t222 - t292;
t151 = pkin(4) * t194 - qJ(5) * t193;
t275 = -t151 + t277;
t274 = t230 * t222 + t249 * t284 + t281;
t197 = (-pkin(4) * t242 + qJ(5) * t243) * t248;
t273 = -t197 + t276;
t272 = rSges(3,1) * t248 - rSges(3,2) * t245;
t271 = -rSges(4,2) * t248 + rSges(4,3) * t245;
t228 = pkin(1) * t249 + t246 * pkin(7);
t263 = -V_base(4) * pkin(6) + t237 * t228 + V_base(2);
t262 = V_base(4) * t227 - t228 * V_base(5) + V_base(3);
t261 = t230 * t292 + t249 * t283 + t274;
t258 = t237 * t199 + t246 * t284 + t263;
t257 = qJD(5) * t191 + t230 * t197 + t261;
t256 = -qJD(3) * t248 + t231 * t198 + t262;
t255 = t237 * t204 + t246 * t283 + t258;
t254 = qJD(4) * t245 + t231 * t205 + t256;
t253 = -qJD(5) * t193 + t237 * t150 + t255;
t252 = qJD(5) * t248 * t243 + t231 * t151 + t254;
t247 = cos(qJ(6));
t244 = sin(qJ(6));
t240 = Icges(2,4) * t249;
t226 = rSges(2,1) * t249 - t246 * rSges(2,2);
t225 = t246 * rSges(2,1) + rSges(2,2) * t249;
t224 = rSges(3,1) * t245 + rSges(3,2) * t248;
t223 = -rSges(4,2) * t245 - rSges(4,3) * t248;
t221 = -qJD(6) * t245 + t237;
t220 = Icges(2,1) * t249 - t297;
t219 = Icges(2,1) * t246 + t240;
t217 = -Icges(2,2) * t246 + t240;
t216 = Icges(2,2) * t249 + t297;
t208 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t207 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t206 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t203 = -pkin(5) * t242 * t248 - pkin(8) * t245;
t196 = -t249 * t282 + t231;
t195 = -t246 * t282 + t230;
t186 = (-t242 * t247 + t243 * t244) * t248;
t185 = (t242 * t244 + t243 * t247) * t248;
t184 = -rSges(4,1) * t249 + t246 * t271;
t183 = t246 * rSges(4,1) + t249 * t271;
t182 = t246 * rSges(3,3) + t249 * t272;
t181 = -rSges(3,3) * t249 + t246 * t272;
t168 = rSges(5,3) * t245 + (-rSges(5,1) * t242 - rSges(5,2) * t243) * t248;
t167 = rSges(6,2) * t245 + (-rSges(6,1) * t242 + rSges(6,3) * t243) * t248;
t156 = pkin(5) * t194 - pkin(8) * t288;
t155 = t192 * pkin(5) - pkin(8) * t287;
t154 = V_base(5) * rSges(2,3) - t225 * t237 + t281;
t153 = t226 * t237 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t152 = t225 * V_base(4) - t226 * V_base(5) + V_base(3);
t149 = -t193 * t244 + t194 * t247;
t148 = -t193 * t247 - t194 * t244;
t147 = t191 * t244 + t192 * t247;
t146 = t191 * t247 - t192 * t244;
t144 = rSges(5,1) * t194 + rSges(5,2) * t193 + rSges(5,3) * t288;
t143 = rSges(6,1) * t194 + rSges(6,2) * t288 - rSges(6,3) * t193;
t142 = t192 * rSges(5,1) - t191 * rSges(5,2) + rSges(5,3) * t287;
t141 = t192 * rSges(6,1) + rSges(6,2) * t287 + t191 * rSges(6,3);
t127 = rSges(7,1) * t186 + rSges(7,2) * t185 - rSges(7,3) * t245;
t126 = Icges(7,1) * t186 + Icges(7,4) * t185 - Icges(7,5) * t245;
t125 = Icges(7,4) * t186 + Icges(7,2) * t185 - Icges(7,6) * t245;
t124 = Icges(7,5) * t186 + Icges(7,6) * t185 - Icges(7,3) * t245;
t123 = t224 * t230 + (-t181 - t227) * t237 + t281;
t122 = t182 * t237 - t224 * t231 + t263;
t121 = t181 * t231 - t182 * t230 + t262;
t120 = rSges(7,1) * t149 + rSges(7,2) * t148 - rSges(7,3) * t288;
t119 = t147 * rSges(7,1) + t146 * rSges(7,2) - rSges(7,3) * t287;
t118 = Icges(7,1) * t149 + Icges(7,4) * t148 - Icges(7,5) * t288;
t117 = Icges(7,1) * t147 + Icges(7,4) * t146 - Icges(7,5) * t287;
t116 = Icges(7,4) * t149 + Icges(7,2) * t148 - Icges(7,6) * t288;
t115 = Icges(7,4) * t147 + Icges(7,2) * t146 - Icges(7,6) * t287;
t114 = Icges(7,5) * t149 + Icges(7,6) * t148 - Icges(7,3) * t288;
t113 = Icges(7,5) * t147 + Icges(7,6) * t146 - Icges(7,3) * t287;
t112 = t223 * t230 + (-t184 + t286) * t237 + t274;
t111 = t183 * t237 + (-t222 - t223) * t231 + t258;
t110 = t184 * t231 + (-t183 - t199) * t230 + t256;
t109 = t168 * t230 + (-t144 + t277) * t237 + t261;
t108 = t142 * t237 + (-t168 + t276) * t231 + t255;
t107 = t144 * t231 + (-t142 + t285) * t230 + t254;
t106 = t167 * t230 + (-t143 + t275) * t237 + t257;
t105 = t141 * t237 + (-t167 + t273) * t231 + t253;
t104 = t143 * t231 + (-t141 + t278) * t230 + t252;
t103 = -t120 * t221 + t127 * t195 + t203 * t230 + (-t156 + t275) * t237 + t257;
t102 = t119 * t221 - t127 * t196 + t155 * t237 + (-t203 + t273) * t231 + t253;
t101 = -t119 * t195 + t120 * t196 + t156 * t231 + (-t155 + t278) * t230 + t252;
t1 = t196 * ((-t113 * t287 + t146 * t115 + t147 * t117) * t196 + (-t114 * t287 + t146 * t116 + t147 * t118) * t195 + (-t124 * t287 + t146 * t125 + t147 * t126) * t221) / 0.2e1 + t195 * ((-t113 * t288 + t115 * t148 + t117 * t149) * t196 + (-t114 * t288 + t148 * t116 + t149 * t118) * t195 + (-t124 * t288 + t125 * t148 + t126 * t149) * t221) / 0.2e1 + t221 * ((-t113 * t245 + t115 * t185 + t117 * t186) * t196 + (-t114 * t245 + t116 * t185 + t118 * t186) * t195 + (-t245 * t124 + t185 * t125 + t186 * t126) * t221) / 0.2e1 + m(1) * (t206 ^ 2 + t207 ^ 2 + t208 ^ 2) / 0.2e1 + m(2) * (t152 ^ 2 + t153 ^ 2 + t154 ^ 2) / 0.2e1 + m(3) * (t121 ^ 2 + t122 ^ 2 + t123 ^ 2) / 0.2e1 + m(4) * (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + m(7) * (t101 ^ 2 + t102 ^ 2 + t103 ^ 2) / 0.2e1 + m(6) * (t104 ^ 2 + t105 ^ 2 + t106 ^ 2) / 0.2e1 + m(5) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + ((-t246 * t216 + t219 * t249 + Icges(1,4)) * V_base(5) + (-t246 * t217 + t249 * t220 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t249 * t216 + t246 * t219 + Icges(1,2)) * V_base(5) + (t217 * t249 + t246 * t220 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (-t301 * t249 + t302 * t246 + (-t193 * t305 + t194 * t303 + t288 * t304) * t237 + (-t193 * t311 + t307 * t194 + t309 * t288) * t231 + (-t310 * t193 + t306 * t194 + t308 * t288) * t230) * t230 / 0.2e1 + (t302 * t249 + t301 * t246 + (t191 * t305 + t192 * t303 + t287 * t304) * t237 + (t311 * t191 + t307 * t192 + t309 * t287) * t231 + (t191 * t310 + t192 * t306 + t308 * t287) * t230) * t231 / 0.2e1 + (((-t307 * t242 + t243 * t311 - t317) * t248 + (t309 + t316) * t245) * t231 + ((-t242 * t306 + t243 * t310 - t318) * t248 + (t308 + t315) * t245) * t230 + (Icges(2,3) + (-t242 * t303 + t243 * t305 - t314) * t248 + (t304 + t313) * t245) * t237) * t237 / 0.2e1 + t237 * V_base(4) * (Icges(2,5) * t249 - Icges(2,6) * t246) + V_base(5) * t237 * (Icges(2,5) * t246 + Icges(2,6) * t249) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
