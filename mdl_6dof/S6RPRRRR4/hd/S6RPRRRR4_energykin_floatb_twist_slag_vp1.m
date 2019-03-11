% Calculate kinetic energy for
% S6RPRRRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR4_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRRR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR4_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR4_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRR4_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:04:43
% EndTime: 2019-03-09 07:04:45
% DurationCPUTime: 2.55s
% Computational Cost: add. (2191->341), mult. (1790->503), div. (0->0), fcn. (1572->12), ass. (0->183)
t245 = sin(pkin(11));
t320 = pkin(2) * t245;
t244 = pkin(11) + qJ(3);
t231 = sin(t244);
t319 = pkin(3) * t231;
t233 = qJ(4) + t244;
t227 = sin(t233);
t318 = pkin(4) * t227;
t246 = cos(pkin(11));
t317 = t246 * pkin(2);
t249 = sin(qJ(1));
t316 = Icges(2,4) * t249;
t315 = Icges(3,4) * t245;
t314 = Icges(3,4) * t246;
t313 = Icges(4,4) * t231;
t232 = cos(t244);
t312 = Icges(4,4) * t232;
t311 = Icges(5,4) * t227;
t228 = cos(t233);
t310 = Icges(5,4) * t228;
t230 = qJ(5) + t233;
t223 = sin(t230);
t309 = Icges(6,4) * t223;
t224 = cos(t230);
t308 = Icges(6,4) * t224;
t307 = t223 * t249;
t251 = cos(qJ(1));
t306 = t223 * t251;
t248 = sin(qJ(6));
t305 = t248 * t251;
t304 = t249 * t248;
t250 = cos(qJ(6));
t303 = t249 * t250;
t302 = t250 * t251;
t156 = -pkin(7) * t251 + t249 * t317;
t217 = t249 * pkin(1) - qJ(2) * t251;
t300 = -t156 - t217;
t299 = pkin(4) * t228;
t298 = pkin(3) * t232;
t295 = qJD(6) * t223;
t294 = -qJD(3) - qJD(4);
t293 = V_base(4) * t217 + V_base(3);
t292 = V_base(5) * pkin(6) + V_base(1);
t127 = -pkin(8) * t251 + t249 * t298;
t289 = -t127 + t300;
t222 = qJD(3) * t249 + V_base(4);
t234 = V_base(6) + qJD(1);
t288 = qJD(2) * t249 + t292;
t123 = -pkin(9) * t251 + t249 * t299;
t287 = -t123 + t289;
t201 = qJD(4) * t249 + t222;
t286 = V_base(5) * t320 + t288;
t285 = pkin(5) * t224 + pkin(10) * t223;
t284 = rSges(3,1) * t246 - rSges(3,2) * t245;
t283 = rSges(4,1) * t232 - rSges(4,2) * t231;
t282 = rSges(5,1) * t228 - rSges(5,2) * t227;
t281 = rSges(6,1) * t224 - rSges(6,2) * t223;
t188 = qJD(5) * t249 + t201;
t280 = Icges(3,1) * t246 - t315;
t279 = Icges(4,1) * t232 - t313;
t278 = Icges(5,1) * t228 - t311;
t277 = Icges(6,1) * t224 - t309;
t276 = -Icges(3,2) * t245 + t314;
t275 = -Icges(4,2) * t231 + t312;
t274 = -Icges(5,2) * t227 + t310;
t273 = -Icges(6,2) * t223 + t308;
t272 = Icges(3,5) * t246 - Icges(3,6) * t245;
t271 = Icges(4,5) * t232 - Icges(4,6) * t231;
t270 = Icges(5,5) * t228 - Icges(5,6) * t227;
t269 = Icges(6,5) * t224 - Icges(6,6) * t223;
t219 = pkin(1) * t251 + t249 * qJ(2);
t268 = -qJD(2) * t251 + t234 * t219 + V_base(2);
t221 = -qJD(3) * t251 + V_base(5);
t267 = t221 * t319 + t286;
t200 = t251 * t294 + V_base(5);
t266 = t200 * t318 + t267;
t187 = V_base(5) + (-qJD(5) + t294) * t251;
t265 = (-Icges(6,3) * t251 + t249 * t269) * t187 + (Icges(6,3) * t249 + t251 * t269) * t188 + (Icges(6,5) * t223 + Icges(6,6) * t224) * t234;
t264 = (-Icges(5,3) * t251 + t249 * t270) * t200 + (Icges(5,3) * t249 + t251 * t270) * t201 + (Icges(5,5) * t227 + Icges(5,6) * t228) * t234;
t263 = (-Icges(4,3) * t251 + t249 * t271) * t221 + (Icges(4,3) * t249 + t251 * t271) * t222 + (Icges(4,5) * t231 + Icges(4,6) * t232) * t234;
t157 = pkin(7) * t249 + t251 * t317;
t262 = V_base(4) * t156 + (-t157 - t219) * V_base(5) + t293;
t261 = (-Icges(3,3) * t251 + t249 * t272) * V_base(5) + (Icges(3,3) * t249 + t251 * t272) * V_base(4) + (Icges(3,5) * t245 + Icges(3,6) * t246) * t234;
t128 = pkin(8) * t249 + t251 * t298;
t260 = t222 * t127 - t128 * t221 + t262;
t259 = t234 * t157 + (-pkin(6) - t320) * V_base(4) + t268;
t124 = pkin(9) * t249 + t251 * t299;
t258 = t201 * t123 - t124 * t200 + t260;
t257 = t234 * t128 - t222 * t319 + t259;
t256 = t234 * t124 - t201 * t318 + t257;
t138 = -Icges(6,6) * t251 + t249 * t273;
t139 = Icges(6,6) * t249 + t251 * t273;
t140 = -Icges(6,5) * t251 + t249 * t277;
t141 = Icges(6,5) * t249 + t251 * t277;
t183 = Icges(6,2) * t224 + t309;
t184 = Icges(6,1) * t223 + t308;
t255 = (-t139 * t223 + t141 * t224) * t188 + (-t138 * t223 + t140 * t224) * t187 + (-t183 * t223 + t184 * t224) * t234;
t148 = -Icges(5,6) * t251 + t249 * t274;
t149 = Icges(5,6) * t249 + t251 * t274;
t150 = -Icges(5,5) * t251 + t249 * t278;
t151 = Icges(5,5) * t249 + t251 * t278;
t190 = Icges(5,2) * t228 + t311;
t191 = Icges(5,1) * t227 + t310;
t254 = (-t149 * t227 + t151 * t228) * t201 + (-t148 * t227 + t150 * t228) * t200 + (-t190 * t227 + t191 * t228) * t234;
t160 = -Icges(4,6) * t251 + t249 * t275;
t161 = Icges(4,6) * t249 + t251 * t275;
t162 = -Icges(4,5) * t251 + t249 * t279;
t163 = Icges(4,5) * t249 + t251 * t279;
t196 = Icges(4,2) * t232 + t313;
t197 = Icges(4,1) * t231 + t312;
t253 = (-t161 * t231 + t163 * t232) * t222 + (-t160 * t231 + t162 * t232) * t221 + (-t196 * t231 + t197 * t232) * t234;
t172 = -Icges(3,6) * t251 + t249 * t276;
t173 = Icges(3,6) * t249 + t251 * t276;
t174 = -Icges(3,5) * t251 + t249 * t280;
t175 = Icges(3,5) * t249 + t251 * t280;
t208 = Icges(3,2) * t246 + t315;
t209 = Icges(3,1) * t245 + t314;
t252 = (-t173 * t245 + t175 * t246) * V_base(4) + (-t172 * t245 + t174 * t246) * V_base(5) + (-t208 * t245 + t209 * t246) * t234;
t241 = Icges(2,4) * t251;
t220 = rSges(2,1) * t251 - t249 * rSges(2,2);
t218 = t249 * rSges(2,1) + rSges(2,2) * t251;
t216 = Icges(2,1) * t251 - t316;
t215 = Icges(2,1) * t249 + t241;
t214 = -Icges(2,2) * t249 + t241;
t213 = Icges(2,2) * t251 + t316;
t212 = Icges(2,5) * t251 - Icges(2,6) * t249;
t211 = Icges(2,5) * t249 + Icges(2,6) * t251;
t210 = rSges(3,1) * t245 + rSges(3,2) * t246;
t206 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t205 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t204 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t199 = rSges(4,1) * t231 + rSges(4,2) * t232;
t198 = -qJD(6) * t224 + t234;
t192 = rSges(5,1) * t227 + rSges(5,2) * t228;
t186 = pkin(5) * t223 - pkin(10) * t224;
t185 = rSges(6,1) * t223 + rSges(6,2) * t224;
t181 = t249 * rSges(3,3) + t251 * t284;
t180 = -rSges(3,3) * t251 + t249 * t284;
t179 = t224 * t302 + t304;
t178 = -t224 * t305 + t303;
t177 = t224 * t303 - t305;
t176 = -t224 * t304 - t302;
t167 = t249 * rSges(4,3) + t251 * t283;
t166 = -rSges(4,3) * t251 + t249 * t283;
t165 = t285 * t251;
t164 = t285 * t249;
t155 = t249 * rSges(5,3) + t251 * t282;
t154 = -rSges(5,3) * t251 + t249 * t282;
t153 = V_base(5) * rSges(2,3) - t218 * t234 + t292;
t152 = t220 * t234 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t145 = t218 * V_base(4) - t220 * V_base(5) + V_base(3);
t143 = t249 * rSges(6,3) + t251 * t281;
t142 = -rSges(6,3) * t251 + t249 * t281;
t135 = t251 * t295 + t188;
t134 = t249 * t295 + t187;
t132 = -rSges(7,3) * t224 + (rSges(7,1) * t250 - rSges(7,2) * t248) * t223;
t131 = -Icges(7,5) * t224 + (Icges(7,1) * t250 - Icges(7,4) * t248) * t223;
t130 = -Icges(7,6) * t224 + (Icges(7,4) * t250 - Icges(7,2) * t248) * t223;
t129 = -Icges(7,3) * t224 + (Icges(7,5) * t250 - Icges(7,6) * t248) * t223;
t122 = t179 * rSges(7,1) + t178 * rSges(7,2) + rSges(7,3) * t306;
t121 = rSges(7,1) * t177 + rSges(7,2) * t176 + rSges(7,3) * t307;
t120 = Icges(7,1) * t179 + Icges(7,4) * t178 + Icges(7,5) * t306;
t119 = Icges(7,1) * t177 + Icges(7,4) * t176 + Icges(7,5) * t307;
t118 = Icges(7,4) * t179 + Icges(7,2) * t178 + Icges(7,6) * t306;
t117 = Icges(7,4) * t177 + Icges(7,2) * t176 + Icges(7,6) * t307;
t116 = Icges(7,5) * t179 + Icges(7,6) * t178 + Icges(7,3) * t306;
t115 = Icges(7,5) * t177 + Icges(7,6) * t176 + Icges(7,3) * t307;
t113 = t210 * V_base(5) + (-t180 - t217) * t234 + t288;
t112 = t234 * t181 + (-pkin(6) - t210) * V_base(4) + t268;
t110 = t180 * V_base(4) + (-t181 - t219) * V_base(5) + t293;
t109 = t199 * t221 + (-t166 + t300) * t234 + t286;
t108 = t234 * t167 - t222 * t199 + t259;
t107 = t166 * t222 - t167 * t221 + t262;
t106 = t192 * t200 + (-t154 + t289) * t234 + t267;
t105 = t234 * t155 - t201 * t192 + t257;
t104 = t154 * t201 - t155 * t200 + t260;
t103 = t185 * t187 + (-t142 + t287) * t234 + t266;
t102 = t234 * t143 - t188 * t185 + t256;
t101 = -t121 * t198 + t132 * t134 + t186 * t187 + (-t164 + t287) * t234 + t266;
t100 = t198 * t122 - t135 * t132 + t234 * t165 - t188 * t186 + t256;
t99 = t142 * t188 - t143 * t187 + t258;
t98 = t121 * t135 - t122 * t134 + t164 * t188 - t165 * t187 + t258;
t1 = t135 * ((t116 * t306 + t178 * t118 + t179 * t120) * t135 + (t115 * t306 + t178 * t117 + t179 * t119) * t134 + (t129 * t306 + t178 * t130 + t179 * t131) * t198) / 0.2e1 + t134 * ((t116 * t307 + t118 * t176 + t120 * t177) * t135 + (t115 * t307 + t176 * t117 + t177 * t119) * t134 + (t129 * t307 + t130 * t176 + t131 * t177) * t198) / 0.2e1 + t201 * (t264 * t249 + t254 * t251) / 0.2e1 + t200 * (t254 * t249 - t264 * t251) / 0.2e1 + t188 * (t265 * t249 + t255 * t251) / 0.2e1 + t187 * (t255 * t249 - t265 * t251) / 0.2e1 + t222 * (t263 * t249 + t253 * t251) / 0.2e1 + t221 * (t253 * t249 - t263 * t251) / 0.2e1 + m(2) * (t145 ^ 2 + t152 ^ 2 + t153 ^ 2) / 0.2e1 + m(4) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + m(3) * (t110 ^ 2 + t112 ^ 2 + t113 ^ 2) / 0.2e1 + m(7) * (t100 ^ 2 + t101 ^ 2 + t98 ^ 2) / 0.2e1 + m(6) * (t102 ^ 2 + t103 ^ 2 + t99 ^ 2) / 0.2e1 + m(5) * (t104 ^ 2 + t105 ^ 2 + t106 ^ 2) / 0.2e1 + m(1) * (t204 ^ 2 + t205 ^ 2 + t206 ^ 2) / 0.2e1 + t198 * ((-t115 * t134 - t116 * t135 - t129 * t198) * t224 + ((-t118 * t248 + t120 * t250) * t135 + (-t117 * t248 + t119 * t250) * t134 + (-t130 * t248 + t131 * t250) * t198) * t223) / 0.2e1 + (t212 * t234 + t261 * t249 + t252 * t251 + (-t249 * t213 + t215 * t251 + Icges(1,4)) * V_base(5) + (-t249 * t214 + t251 * t216 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t211 * t234 + t252 * t249 - t261 * t251 + (t251 * t213 + t249 * t215 + Icges(1,2)) * V_base(5) + (t214 * t251 + t249 * t216 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t161 * t232 + t163 * t231) * t222 + (t160 * t232 + t162 * t231) * t221 + (t149 * t228 + t151 * t227) * t201 + (t148 * t228 + t150 * t227) * t200 + (t139 * t224 + t141 * t223) * t188 + (t138 * t224 + t140 * t223) * t187 + (t172 * t246 + t174 * t245 + t211) * V_base(5) + (t173 * t246 + t175 * t245 + t212) * V_base(4) + (t224 * t183 + t223 * t184 + t228 * t190 + t227 * t191 + t232 * t196 + t231 * t197 + t246 * t208 + t245 * t209 + Icges(2,3)) * t234) * t234 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
