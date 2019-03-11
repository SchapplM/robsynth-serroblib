% Calculate kinetic energy for
% S6RRRRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR4_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRPR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR4_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR4_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPR4_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:06:06
% EndTime: 2019-03-09 22:06:10
% DurationCPUTime: 3.86s
% Computational Cost: add. (2463->380), mult. (2382->561), div. (0->0), fcn. (2276->12), ass. (0->181)
t332 = -Icges(6,3) - Icges(5,3);
t260 = qJ(4) + pkin(11);
t248 = sin(t260);
t249 = cos(t260);
t268 = cos(qJ(1));
t261 = qJ(2) + qJ(3);
t255 = cos(t261);
t265 = sin(qJ(1));
t311 = t255 * t265;
t189 = -t248 * t311 - t249 * t268;
t190 = -t248 * t268 + t249 * t311;
t266 = cos(qJ(4));
t306 = t266 * t268;
t263 = sin(qJ(4));
t309 = t263 * t265;
t204 = -t255 * t309 - t306;
t307 = t265 * t266;
t308 = t263 * t268;
t205 = t255 * t307 - t308;
t254 = sin(t261);
t313 = t254 * t265;
t331 = Icges(5,5) * t205 + Icges(6,5) * t190 + Icges(5,6) * t204 + Icges(6,6) * t189 - t313 * t332;
t310 = t255 * t268;
t191 = -t248 * t310 + t249 * t265;
t192 = t248 * t265 + t249 * t310;
t206 = -t255 * t308 + t307;
t207 = t255 * t306 + t309;
t312 = t254 * t268;
t330 = Icges(5,5) * t207 + Icges(6,5) * t192 + Icges(5,6) * t206 + Icges(6,6) * t191 - t312 * t332;
t329 = t332 * t255 + (Icges(5,5) * t266 + Icges(6,5) * t249 - Icges(5,6) * t263 - Icges(6,6) * t248) * t254;
t264 = sin(qJ(2));
t324 = pkin(2) * t264;
t267 = cos(qJ(2));
t322 = pkin(2) * t267;
t321 = t266 * pkin(4);
t318 = Icges(2,4) * t265;
t317 = Icges(3,4) * t264;
t316 = Icges(3,4) * t267;
t315 = Icges(4,4) * t254;
t314 = Icges(4,4) * t255;
t176 = -pkin(8) * t268 + t265 * t322;
t240 = t265 * pkin(1) - t268 * pkin(7);
t305 = -t176 - t240;
t304 = pkin(5) * t249;
t302 = qJD(4) * t254;
t301 = qJD(5) * t254;
t300 = qJD(6) * t254;
t299 = V_base(5) * pkin(6) + V_base(1);
t243 = qJD(2) * t265 + V_base(4);
t251 = V_base(6) + qJD(1);
t296 = pkin(5) * t248;
t242 = -qJD(2) * t268 + V_base(5);
t295 = t242 * t324 + t299;
t215 = qJD(3) * t265 + t243;
t294 = pkin(3) * t255 + pkin(9) * t254;
t293 = rSges(3,1) * t267 - rSges(3,2) * t264;
t292 = rSges(4,1) * t255 - rSges(4,2) * t254;
t188 = t268 * t302 + t215;
t291 = Icges(3,1) * t267 - t317;
t290 = Icges(4,1) * t255 - t315;
t289 = -Icges(3,2) * t264 + t316;
t288 = -Icges(4,2) * t254 + t314;
t287 = Icges(3,5) * t267 - Icges(3,6) * t264;
t286 = Icges(4,5) * t255 - Icges(4,6) * t254;
t241 = t268 * pkin(1) + t265 * pkin(7);
t285 = -V_base(4) * pkin(6) + t251 * t241 + V_base(2);
t284 = V_base(4) * t240 - t241 * V_base(5) + V_base(3);
t214 = V_base(5) + (-qJD(2) - qJD(3)) * t268;
t187 = t265 * t302 + t214;
t283 = qJ(5) * t254 + t255 * t321;
t282 = (-Icges(4,3) * t268 + t265 * t286) * t214 + (Icges(4,3) * t265 + t268 * t286) * t215 + (Icges(4,5) * t254 + Icges(4,6) * t255) * t251;
t281 = (-Icges(3,3) * t268 + t265 * t287) * t242 + (Icges(3,3) * t265 + t268 * t287) * t243 + (Icges(3,5) * t264 + Icges(3,6) * t267) * t251;
t280 = pkin(10) * t254 + t255 * t304;
t201 = t294 * t265;
t213 = pkin(3) * t254 - pkin(9) * t255;
t279 = t214 * t213 + (-t201 + t305) * t251 + t295;
t177 = pkin(8) * t265 + t268 * t322;
t278 = t243 * t176 - t177 * t242 + t284;
t277 = t251 * t177 - t243 * t324 + t285;
t154 = -qJ(5) * t255 + t254 * t321;
t276 = t187 * t154 + t268 * t301 + t279;
t202 = t294 * t268;
t275 = t215 * t201 - t202 * t214 + t278;
t274 = t251 * t202 - t213 * t215 + t277;
t140 = pkin(4) * t309 + t268 * t283;
t222 = -qJD(4) * t255 + t251;
t273 = t222 * t140 + t265 * t301 + t274;
t139 = -pkin(4) * t308 + t265 * t283;
t272 = -qJD(5) * t255 + t188 * t139 + t275;
t181 = -Icges(4,6) * t268 + t265 * t288;
t182 = Icges(4,6) * t265 + t268 * t288;
t183 = -Icges(4,5) * t268 + t265 * t290;
t184 = Icges(4,5) * t265 + t268 * t290;
t210 = Icges(4,2) * t255 + t315;
t211 = Icges(4,1) * t254 + t314;
t271 = (-t182 * t254 + t184 * t255) * t215 + (-t181 * t254 + t183 * t255) * t214 + (-t210 * t254 + t211 * t255) * t251;
t195 = -Icges(3,6) * t268 + t265 * t289;
t196 = Icges(3,6) * t265 + t268 * t289;
t197 = -Icges(3,5) * t268 + t265 * t291;
t198 = Icges(3,5) * t265 + t268 * t291;
t227 = Icges(3,2) * t267 + t317;
t230 = Icges(3,1) * t264 + t316;
t270 = (-t196 * t264 + t198 * t267) * t243 + (-t195 * t264 + t197 * t267) * t242 + (-t227 * t264 + t230 * t267) * t251;
t256 = Icges(2,4) * t268;
t250 = qJ(6) + t260;
t245 = cos(t250);
t244 = sin(t250);
t235 = rSges(2,1) * t268 - rSges(2,2) * t265;
t234 = rSges(2,1) * t265 + rSges(2,2) * t268;
t233 = rSges(3,1) * t264 + rSges(3,2) * t267;
t232 = Icges(2,1) * t268 - t318;
t231 = Icges(2,1) * t265 + t256;
t229 = -Icges(2,2) * t265 + t256;
t228 = Icges(2,2) * t268 + t318;
t221 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t220 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t219 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t212 = rSges(4,1) * t254 + rSges(4,2) * t255;
t203 = (-qJD(4) - qJD(6)) * t255 + t251;
t200 = rSges(3,3) * t265 + t268 * t293;
t199 = -rSges(3,3) * t268 + t265 * t293;
t186 = rSges(4,3) * t265 + t268 * t292;
t185 = -rSges(4,3) * t268 + t265 * t292;
t175 = t244 * t265 + t245 * t310;
t174 = -t244 * t310 + t245 * t265;
t173 = -t244 * t268 + t245 * t311;
t172 = -t244 * t311 - t245 * t268;
t171 = -rSges(5,3) * t255 + (rSges(5,1) * t266 - rSges(5,2) * t263) * t254;
t170 = -Icges(5,5) * t255 + (Icges(5,1) * t266 - Icges(5,4) * t263) * t254;
t169 = -Icges(5,6) * t255 + (Icges(5,4) * t266 - Icges(5,2) * t263) * t254;
t167 = V_base(5) * rSges(2,3) - t234 * t251 + t299;
t166 = t235 * t251 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t164 = t234 * V_base(4) - t235 * V_base(5) + V_base(3);
t162 = -rSges(6,3) * t255 + (rSges(6,1) * t249 - rSges(6,2) * t248) * t254;
t161 = -Icges(6,5) * t255 + (Icges(6,1) * t249 - Icges(6,4) * t248) * t254;
t160 = -Icges(6,6) * t255 + (Icges(6,4) * t249 - Icges(6,2) * t248) * t254;
t158 = t268 * t300 + t188;
t157 = t265 * t300 + t187;
t155 = -rSges(7,3) * t255 + (rSges(7,1) * t245 - rSges(7,2) * t244) * t254;
t153 = -Icges(7,5) * t255 + (Icges(7,1) * t245 - Icges(7,4) * t244) * t254;
t152 = -Icges(7,6) * t255 + (Icges(7,4) * t245 - Icges(7,2) * t244) * t254;
t151 = -Icges(7,3) * t255 + (Icges(7,5) * t245 - Icges(7,6) * t244) * t254;
t149 = -pkin(10) * t255 + t254 * t304;
t148 = rSges(5,1) * t207 + rSges(5,2) * t206 + rSges(5,3) * t312;
t147 = rSges(5,1) * t205 + rSges(5,2) * t204 + rSges(5,3) * t313;
t146 = Icges(5,1) * t207 + Icges(5,4) * t206 + Icges(5,5) * t312;
t145 = Icges(5,1) * t205 + Icges(5,4) * t204 + Icges(5,5) * t313;
t144 = Icges(5,4) * t207 + Icges(5,2) * t206 + Icges(5,6) * t312;
t143 = Icges(5,4) * t205 + Icges(5,2) * t204 + Icges(5,6) * t313;
t137 = rSges(6,1) * t192 + rSges(6,2) * t191 + rSges(6,3) * t312;
t136 = rSges(6,1) * t190 + rSges(6,2) * t189 + rSges(6,3) * t313;
t135 = Icges(6,1) * t192 + Icges(6,4) * t191 + Icges(6,5) * t312;
t134 = Icges(6,1) * t190 + Icges(6,4) * t189 + Icges(6,5) * t313;
t133 = Icges(6,4) * t192 + Icges(6,2) * t191 + Icges(6,6) * t312;
t132 = Icges(6,4) * t190 + Icges(6,2) * t189 + Icges(6,6) * t313;
t128 = rSges(7,1) * t175 + rSges(7,2) * t174 + rSges(7,3) * t312;
t127 = rSges(7,1) * t173 + rSges(7,2) * t172 + rSges(7,3) * t313;
t126 = Icges(7,1) * t175 + Icges(7,4) * t174 + Icges(7,5) * t312;
t125 = Icges(7,1) * t173 + Icges(7,4) * t172 + Icges(7,5) * t313;
t124 = Icges(7,4) * t175 + Icges(7,2) * t174 + Icges(7,6) * t312;
t123 = Icges(7,4) * t173 + Icges(7,2) * t172 + Icges(7,6) * t313;
t122 = Icges(7,5) * t175 + Icges(7,6) * t174 + Icges(7,3) * t312;
t121 = Icges(7,5) * t173 + Icges(7,6) * t172 + Icges(7,3) * t313;
t120 = t233 * t242 + (-t199 - t240) * t251 + t299;
t119 = t200 * t251 - t233 * t243 + t285;
t117 = t265 * t296 + t268 * t280;
t116 = t265 * t280 - t268 * t296;
t115 = t199 * t243 - t200 * t242 + t284;
t114 = t212 * t214 + (-t185 + t305) * t251 + t295;
t113 = t186 * t251 - t212 * t215 + t277;
t112 = t185 * t215 - t186 * t214 + t278;
t111 = -t147 * t222 + t171 * t187 + t279;
t110 = t148 * t222 - t171 * t188 + t274;
t109 = t147 * t188 - t148 * t187 + t275;
t108 = t162 * t187 + (-t136 - t139) * t222 + t276;
t107 = t137 * t222 + (-t154 - t162) * t188 + t273;
t106 = t136 * t188 + (-t137 - t140) * t187 + t272;
t105 = -t127 * t203 + t149 * t187 + t155 * t157 + (-t116 - t139) * t222 + t276;
t104 = t117 * t222 + t128 * t203 - t155 * t158 + (-t149 - t154) * t188 + t273;
t103 = t116 * t188 + t127 * t158 - t128 * t157 + (-t117 - t140) * t187 + t272;
t1 = t203 * ((-t121 * t157 - t122 * t158 - t151 * t203) * t255 + ((-t124 * t244 + t126 * t245) * t158 + (-t123 * t244 + t125 * t245) * t157 + (-t152 * t244 + t153 * t245) * t203) * t254) / 0.2e1 + m(1) * (t219 ^ 2 + t220 ^ 2 + t221 ^ 2) / 0.2e1 + m(2) * (t164 ^ 2 + t166 ^ 2 + t167 ^ 2) / 0.2e1 + m(3) * (t115 ^ 2 + t119 ^ 2 + t120 ^ 2) / 0.2e1 + m(4) * (t112 ^ 2 + t113 ^ 2 + t114 ^ 2) / 0.2e1 + m(7) * (t103 ^ 2 + t104 ^ 2 + t105 ^ 2) / 0.2e1 + m(6) * (t106 ^ 2 + t107 ^ 2 + t108 ^ 2) / 0.2e1 + m(5) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + t157 * ((t122 * t313 + t124 * t172 + t126 * t173) * t158 + (t121 * t313 + t123 * t172 + t125 * t173) * t157 + (t151 * t313 + t152 * t172 + t153 * t173) * t203) / 0.2e1 + t158 * ((t122 * t312 + t124 * t174 + t126 * t175) * t158 + (t121 * t312 + t123 * t174 + t125 * t175) * t157 + (t151 * t312 + t152 * t174 + t153 * t175) * t203) / 0.2e1 + t243 * (t265 * t281 + t268 * t270) / 0.2e1 + t242 * (t265 * t270 - t268 * t281) / 0.2e1 + t215 * (t265 * t282 + t268 * t271) / 0.2e1 + t214 * (t265 * t271 - t282 * t268) / 0.2e1 + ((t160 * t189 + t161 * t190 + t169 * t204 + t170 * t205 + t313 * t329) * t222 + (t133 * t189 + t135 * t190 + t144 * t204 + t146 * t205 + t313 * t330) * t188 + (t132 * t189 + t134 * t190 + t143 * t204 + t145 * t205 + t331 * t313) * t187) * t187 / 0.2e1 + ((t160 * t191 + t161 * t192 + t169 * t206 + t170 * t207 + t312 * t329) * t222 + (t133 * t191 + t135 * t192 + t144 * t206 + t146 * t207 + t330 * t312) * t188 + (t132 * t191 + t134 * t192 + t143 * t206 + t145 * t207 + t312 * t331) * t187) * t188 / 0.2e1 + ((-t187 * t331 - t330 * t188 - t329 * t222) * t255 + ((-t160 * t248 + t161 * t249 - t169 * t263 + t170 * t266) * t222 + (-t133 * t248 + t135 * t249 - t144 * t263 + t146 * t266) * t188 + (-t132 * t248 + t134 * t249 - t143 * t263 + t145 * t266) * t187) * t254) * t222 / 0.2e1 + ((-t228 * t265 + t231 * t268 + Icges(1,4)) * V_base(5) + (-t229 * t265 + t232 * t268 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t228 * t268 + t231 * t265 + Icges(1,2)) * V_base(5) + (t229 * t268 + t232 * t265 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t196 * t267 + t198 * t264) * t243 + (t195 * t267 + t197 * t264) * t242 + (t182 * t255 + t184 * t254) * t215 + (t181 * t255 + t183 * t254) * t214 + (t210 * t255 + t211 * t254 + t227 * t267 + t230 * t264 + Icges(2,3)) * t251) * t251 / 0.2e1 + t251 * V_base(4) * (Icges(2,5) * t268 - Icges(2,6) * t265) + V_base(5) * t251 * (Icges(2,5) * t265 + Icges(2,6) * t268) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
