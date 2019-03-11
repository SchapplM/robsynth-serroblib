% Calculate kinetic energy for
% S6RRPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
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
% Datum: 2019-03-09 08:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR2_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPPRR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR2_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR2_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRR2_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:49:41
% EndTime: 2019-03-09 08:49:44
% DurationCPUTime: 3.81s
% Computational Cost: add. (2345->385), mult. (2273->554), div. (0->0), fcn. (2167->12), ass. (0->185)
t337 = Icges(3,3) + Icges(4,3);
t261 = qJ(2) + pkin(10);
t249 = sin(t261);
t251 = cos(t261);
t266 = sin(qJ(2));
t268 = cos(qJ(2));
t336 = Icges(3,5) * t268 + Icges(4,5) * t251 - Icges(3,6) * t266 - Icges(4,6) * t249;
t267 = sin(qJ(1));
t269 = cos(qJ(1));
t321 = Icges(4,4) * t251;
t288 = -Icges(4,2) * t249 + t321;
t182 = -Icges(4,6) * t269 + t267 * t288;
t183 = Icges(4,6) * t267 + t269 * t288;
t322 = Icges(4,4) * t249;
t290 = Icges(4,1) * t251 - t322;
t184 = -Icges(4,5) * t269 + t267 * t290;
t185 = Icges(4,5) * t267 + t269 * t290;
t323 = Icges(3,4) * t268;
t289 = -Icges(3,2) * t266 + t323;
t195 = -Icges(3,6) * t269 + t267 * t289;
t196 = Icges(3,6) * t267 + t269 * t289;
t324 = Icges(3,4) * t266;
t291 = Icges(3,1) * t268 - t324;
t197 = -Icges(3,5) * t269 + t267 * t291;
t198 = Icges(3,5) * t267 + t269 * t291;
t212 = Icges(4,2) * t251 + t322;
t213 = Icges(4,1) * t249 + t321;
t227 = Icges(3,2) * t268 + t324;
t230 = Icges(3,1) * t266 + t323;
t242 = -qJD(2) * t269 + V_base(5);
t243 = qJD(2) * t267 + V_base(4);
t253 = V_base(6) + qJD(1);
t335 = (-t212 * t249 + t213 * t251 - t227 * t266 + t230 * t268) * t253 + (-t183 * t249 + t185 * t251 - t196 * t266 + t198 * t268) * t243 + (-t182 * t249 + t184 * t251 - t195 * t266 + t197 * t268) * t242;
t334 = (Icges(3,5) * t266 + Icges(4,5) * t249 + Icges(3,6) * t268 + Icges(4,6) * t251) * t253 + (t267 * t337 + t336 * t269) * t243 + (t336 * t267 - t269 * t337) * t242;
t330 = pkin(2) * t266;
t328 = pkin(2) * t268;
t263 = cos(pkin(11));
t327 = t263 * pkin(4);
t325 = Icges(2,4) * t267;
t320 = t249 * t267;
t319 = t249 * t269;
t318 = t251 * t269;
t262 = sin(pkin(11));
t317 = t262 * t269;
t316 = t263 * t269;
t260 = pkin(11) + qJ(5);
t252 = qJ(6) + t260;
t244 = sin(t252);
t315 = t267 * t244;
t245 = cos(t252);
t314 = t267 * t245;
t248 = sin(t260);
t313 = t267 * t248;
t250 = cos(t260);
t312 = t267 * t250;
t311 = t267 * t262;
t310 = t267 * t263;
t178 = -qJ(3) * t269 + t267 * t328;
t240 = t267 * pkin(1) - pkin(7) * t269;
t308 = -t178 - t240;
t179 = qJ(3) * t267 + t269 * t328;
t292 = pkin(3) * t251 + qJ(4) * t249;
t200 = t292 * t269;
t307 = -t179 - t200;
t306 = pkin(5) * t250;
t304 = qJD(4) * t249;
t303 = qJD(5) * t249;
t302 = qJD(6) * t249;
t301 = V_base(5) * pkin(6) + V_base(1);
t199 = t292 * t267;
t298 = -t199 + t308;
t214 = pkin(3) * t249 - qJ(4) * t251;
t297 = -t214 - t330;
t296 = pkin(5) * t248;
t209 = t269 * t303 + t243;
t295 = qJD(3) * t267 + t242 * t330 + t301;
t294 = rSges(3,1) * t268 - rSges(3,2) * t266;
t293 = rSges(4,1) * t251 - rSges(4,2) * t249;
t208 = t267 * t303 + t242;
t241 = pkin(1) * t269 + t267 * pkin(7);
t285 = -V_base(4) * pkin(6) + t253 * t241 + V_base(2);
t284 = V_base(4) * t240 - t241 * V_base(5) + V_base(3);
t283 = t242 * t214 + t269 * t304 + t295;
t282 = t243 * t178 + t284;
t279 = pkin(8) * t249 + t251 * t327;
t278 = pkin(9) * t249 + t251 * t306;
t277 = -qJD(3) * t269 + t253 * t179 + t285;
t276 = -qJD(4) * t251 + t243 * t199 + t282;
t275 = t253 * t200 + t267 * t304 + t277;
t140 = -pkin(4) * t317 + t267 * t279;
t152 = -pkin(8) * t251 + t249 * t327;
t274 = t242 * t152 + (-t140 + t298) * t253 + t283;
t141 = pkin(4) * t311 + t269 * t279;
t273 = t243 * t140 + (-t141 + t307) * t242 + t276;
t272 = t253 * t141 + (-t152 + t297) * t243 + t275;
t257 = Icges(2,4) * t269;
t239 = rSges(2,1) * t269 - t267 * rSges(2,2);
t238 = t267 * rSges(2,1) + rSges(2,2) * t269;
t237 = rSges(3,1) * t266 + rSges(3,2) * t268;
t232 = Icges(2,1) * t269 - t325;
t231 = Icges(2,1) * t267 + t257;
t229 = -Icges(2,2) * t267 + t257;
t228 = Icges(2,2) * t269 + t325;
t223 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t222 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t221 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t219 = -qJD(5) * t251 + t253;
t215 = rSges(4,1) * t249 + rSges(4,2) * t251;
t207 = t251 * t316 + t311;
t206 = -t251 * t317 + t310;
t205 = t251 * t310 - t317;
t204 = -t251 * t311 - t316;
t203 = (-qJD(5) - qJD(6)) * t251 + t253;
t202 = t267 * rSges(3,3) + t269 * t294;
t201 = -rSges(3,3) * t269 + t267 * t294;
t192 = t250 * t318 + t313;
t191 = -t248 * t318 + t312;
t190 = -t248 * t269 + t251 * t312;
t189 = -t250 * t269 - t251 * t313;
t187 = t267 * rSges(4,3) + t269 * t293;
t186 = -rSges(4,3) * t269 + t267 * t293;
t177 = t245 * t318 + t315;
t176 = -t244 * t318 + t314;
t175 = -t244 * t269 + t251 * t314;
t174 = -t245 * t269 - t251 * t315;
t172 = V_base(5) * rSges(2,3) - t238 * t253 + t301;
t171 = t239 * t253 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t170 = -rSges(5,3) * t251 + (rSges(5,1) * t263 - rSges(5,2) * t262) * t249;
t169 = -Icges(5,5) * t251 + (Icges(5,1) * t263 - Icges(5,4) * t262) * t249;
t168 = -Icges(5,6) * t251 + (Icges(5,4) * t263 - Icges(5,2) * t262) * t249;
t167 = -Icges(5,3) * t251 + (Icges(5,5) * t263 - Icges(5,6) * t262) * t249;
t166 = t269 * t302 + t209;
t165 = t267 * t302 + t208;
t163 = t238 * V_base(4) - t239 * V_base(5) + V_base(3);
t161 = -rSges(6,3) * t251 + (rSges(6,1) * t250 - rSges(6,2) * t248) * t249;
t160 = -Icges(6,5) * t251 + (Icges(6,1) * t250 - Icges(6,4) * t248) * t249;
t159 = -Icges(6,6) * t251 + (Icges(6,4) * t250 - Icges(6,2) * t248) * t249;
t158 = -Icges(6,3) * t251 + (Icges(6,5) * t250 - Icges(6,6) * t248) * t249;
t157 = -rSges(7,3) * t251 + (rSges(7,1) * t245 - rSges(7,2) * t244) * t249;
t155 = -Icges(7,5) * t251 + (Icges(7,1) * t245 - Icges(7,4) * t244) * t249;
t154 = -Icges(7,6) * t251 + (Icges(7,4) * t245 - Icges(7,2) * t244) * t249;
t153 = -Icges(7,3) * t251 + (Icges(7,5) * t245 - Icges(7,6) * t244) * t249;
t150 = -pkin(9) * t251 + t249 * t306;
t149 = t207 * rSges(5,1) + t206 * rSges(5,2) + rSges(5,3) * t319;
t148 = rSges(5,1) * t205 + rSges(5,2) * t204 + rSges(5,3) * t320;
t147 = Icges(5,1) * t207 + Icges(5,4) * t206 + Icges(5,5) * t319;
t146 = Icges(5,1) * t205 + Icges(5,4) * t204 + Icges(5,5) * t320;
t145 = Icges(5,4) * t207 + Icges(5,2) * t206 + Icges(5,6) * t319;
t144 = Icges(5,4) * t205 + Icges(5,2) * t204 + Icges(5,6) * t320;
t143 = Icges(5,5) * t207 + Icges(5,6) * t206 + Icges(5,3) * t319;
t142 = Icges(5,5) * t205 + Icges(5,6) * t204 + Icges(5,3) * t320;
t138 = t192 * rSges(6,1) + t191 * rSges(6,2) + rSges(6,3) * t319;
t137 = rSges(6,1) * t190 + rSges(6,2) * t189 + rSges(6,3) * t320;
t136 = Icges(6,1) * t192 + Icges(6,4) * t191 + Icges(6,5) * t319;
t135 = Icges(6,1) * t190 + Icges(6,4) * t189 + Icges(6,5) * t320;
t134 = Icges(6,4) * t192 + Icges(6,2) * t191 + Icges(6,6) * t319;
t133 = Icges(6,4) * t190 + Icges(6,2) * t189 + Icges(6,6) * t320;
t132 = Icges(6,5) * t192 + Icges(6,6) * t191 + Icges(6,3) * t319;
t131 = Icges(6,5) * t190 + Icges(6,6) * t189 + Icges(6,3) * t320;
t129 = t177 * rSges(7,1) + t176 * rSges(7,2) + rSges(7,3) * t319;
t128 = rSges(7,1) * t175 + rSges(7,2) * t174 + rSges(7,3) * t320;
t127 = Icges(7,1) * t177 + Icges(7,4) * t176 + Icges(7,5) * t319;
t126 = Icges(7,1) * t175 + Icges(7,4) * t174 + Icges(7,5) * t320;
t125 = Icges(7,4) * t177 + Icges(7,2) * t176 + Icges(7,6) * t319;
t124 = Icges(7,4) * t175 + Icges(7,2) * t174 + Icges(7,6) * t320;
t123 = Icges(7,5) * t177 + Icges(7,6) * t176 + Icges(7,3) * t319;
t122 = Icges(7,5) * t175 + Icges(7,6) * t174 + Icges(7,3) * t320;
t121 = t237 * t242 + (-t201 - t240) * t253 + t301;
t120 = t202 * t253 - t237 * t243 + t285;
t119 = t267 * t296 + t269 * t278;
t118 = t267 * t278 - t269 * t296;
t117 = t201 * t243 - t202 * t242 + t284;
t116 = t215 * t242 + (-t186 + t308) * t253 + t295;
t115 = t253 * t187 + (-t215 - t330) * t243 + t277;
t114 = t186 * t243 + (-t179 - t187) * t242 + t282;
t113 = t170 * t242 + (-t148 + t298) * t253 + t283;
t112 = t253 * t149 + (-t170 + t297) * t243 + t275;
t111 = t148 * t243 + (-t149 + t307) * t242 + t276;
t110 = -t137 * t219 + t161 * t208 + t274;
t109 = t219 * t138 - t209 * t161 + t272;
t108 = t137 * t209 - t138 * t208 + t273;
t107 = -t118 * t219 - t128 * t203 + t150 * t208 + t157 * t165 + t274;
t106 = t219 * t119 + t203 * t129 - t209 * t150 - t166 * t157 + t272;
t105 = t118 * t209 - t119 * t208 + t128 * t166 - t129 * t165 + t273;
t1 = t219 * ((-t131 * t208 - t132 * t209 - t158 * t219) * t251 + ((-t134 * t248 + t136 * t250) * t209 + (-t133 * t248 + t135 * t250) * t208 + (-t159 * t248 + t160 * t250) * t219) * t249) / 0.2e1 + t203 * ((-t122 * t165 - t123 * t166 - t153 * t203) * t251 + ((-t125 * t244 + t127 * t245) * t166 + (-t124 * t244 + t126 * t245) * t165 + (-t154 * t244 + t155 * t245) * t203) * t249) / 0.2e1 + m(1) * (t221 ^ 2 + t222 ^ 2 + t223 ^ 2) / 0.2e1 + m(2) * (t163 ^ 2 + t171 ^ 2 + t172 ^ 2) / 0.2e1 + m(3) * (t117 ^ 2 + t120 ^ 2 + t121 ^ 2) / 0.2e1 + m(5) * (t111 ^ 2 + t112 ^ 2 + t113 ^ 2) / 0.2e1 + m(4) * (t114 ^ 2 + t115 ^ 2 + t116 ^ 2) / 0.2e1 + m(7) * (t105 ^ 2 + t106 ^ 2 + t107 ^ 2) / 0.2e1 + m(6) * (t108 ^ 2 + t109 ^ 2 + t110 ^ 2) / 0.2e1 + t209 * ((t132 * t319 + t191 * t134 + t192 * t136) * t209 + (t131 * t319 + t191 * t133 + t192 * t135) * t208 + (t158 * t319 + t191 * t159 + t192 * t160) * t219) / 0.2e1 + t166 * ((t123 * t319 + t176 * t125 + t177 * t127) * t166 + (t122 * t319 + t176 * t124 + t177 * t126) * t165 + (t153 * t319 + t176 * t154 + t177 * t155) * t203) / 0.2e1 + t208 * ((t132 * t320 + t134 * t189 + t136 * t190) * t209 + (t131 * t320 + t133 * t189 + t135 * t190) * t208 + (t158 * t320 + t159 * t189 + t160 * t190) * t219) / 0.2e1 + t165 * ((t123 * t320 + t125 * t174 + t127 * t175) * t166 + (t122 * t320 + t124 * t174 + t126 * t175) * t165 + (t153 * t320 + t154 * t174 + t155 * t175) * t203) / 0.2e1 + ((-t267 * t228 + t231 * t269 + Icges(1,4)) * V_base(5) + (-t267 * t229 + t232 * t269 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t228 * t269 + t267 * t231 + Icges(1,2)) * V_base(5) + (t229 * t269 + t267 * t232 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t143 * t320 + t145 * t204 + t147 * t205) * t243 + (t142 * t320 + t144 * t204 + t146 * t205) * t242 + (t167 * t320 + t168 * t204 + t169 * t205) * t253 - t334 * t269 + t335 * t267) * t242 / 0.2e1 + ((t143 * t319 + t206 * t145 + t207 * t147) * t243 + (t142 * t319 + t206 * t144 + t207 * t146) * t242 + (t167 * t319 + t206 * t168 + t207 * t169) * t253 + t335 * t269 + t334 * t267) * t243 / 0.2e1 + ((t196 * t268 + t198 * t266 + (-t143 + t183) * t251 + (-t145 * t262 + t147 * t263 + t185) * t249) * t243 + (t195 * t268 + t197 * t266 + (-t142 + t182) * t251 + (-t144 * t262 + t146 * t263 + t184) * t249) * t242 + (t227 * t268 + t230 * t266 + Icges(2,3) + (-t167 + t212) * t251 + (-t168 * t262 + t169 * t263 + t213) * t249) * t253) * t253 / 0.2e1 + t253 * V_base(4) * (Icges(2,5) * t269 - Icges(2,6) * t267) + t253 * V_base(5) * (Icges(2,5) * t267 + Icges(2,6) * t269) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
