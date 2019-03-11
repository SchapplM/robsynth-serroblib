% Calculate kinetic energy for
% S6RRRRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 03:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR2_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR2_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRRR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR2_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR2_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR2_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRR2_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:33:06
% EndTime: 2019-03-10 03:33:09
% DurationCPUTime: 2.84s
% Computational Cost: add. (2424->362), mult. (2114->551), div. (0->0), fcn. (1952->12), ass. (0->187)
t265 = sin(qJ(2));
t335 = pkin(2) * t265;
t263 = qJ(2) + qJ(3);
t254 = sin(t263);
t334 = pkin(3) * t254;
t268 = cos(qJ(2));
t333 = t268 * pkin(2);
t267 = cos(qJ(5));
t332 = pkin(5) * t267;
t266 = sin(qJ(1));
t329 = Icges(2,4) * t266;
t328 = Icges(3,4) * t265;
t327 = Icges(3,4) * t268;
t326 = Icges(4,4) * t254;
t256 = cos(t263);
t325 = Icges(4,4) * t256;
t258 = qJ(4) + t263;
t245 = sin(t258);
t324 = Icges(5,4) * t245;
t246 = cos(t258);
t323 = Icges(5,4) * t246;
t322 = t245 * t266;
t269 = cos(qJ(1));
t321 = t245 * t269;
t262 = qJ(5) + qJ(6);
t253 = sin(t262);
t320 = t253 * t266;
t319 = t253 * t269;
t255 = cos(t262);
t318 = t255 * t266;
t317 = t255 * t269;
t264 = sin(qJ(5));
t316 = t264 * t266;
t315 = t264 * t269;
t314 = t266 * t267;
t313 = t267 * t269;
t176 = -pkin(8) * t269 + t266 * t333;
t241 = t266 * pkin(1) - t269 * pkin(7);
t312 = -t176 - t241;
t311 = pkin(3) * t256;
t309 = qJD(5) * t245;
t308 = qJD(6) * t245;
t307 = -qJD(2) - qJD(3);
t306 = V_base(5) * pkin(6) + V_base(1);
t148 = -pkin(9) * t269 + t266 * t311;
t303 = -t148 + t312;
t244 = qJD(2) * t266 + V_base(4);
t249 = V_base(6) + qJD(1);
t243 = -qJD(2) * t269 + V_base(5);
t302 = t243 * t335 + t306;
t219 = qJD(3) * t266 + t244;
t218 = t269 * t307 + V_base(5);
t301 = t218 * t334 + t302;
t300 = pkin(4) * t246 + pkin(10) * t245;
t299 = rSges(3,1) * t268 - rSges(3,2) * t265;
t298 = rSges(4,1) * t256 - rSges(4,2) * t254;
t297 = rSges(5,1) * t246 - rSges(5,2) * t245;
t207 = qJD(4) * t266 + t219;
t296 = Icges(3,1) * t268 - t328;
t295 = Icges(4,1) * t256 - t326;
t294 = Icges(5,1) * t246 - t324;
t293 = -Icges(3,2) * t265 + t327;
t292 = -Icges(4,2) * t254 + t325;
t291 = -Icges(5,2) * t245 + t323;
t290 = Icges(3,5) * t268 - Icges(3,6) * t265;
t289 = Icges(4,5) * t256 - Icges(4,6) * t254;
t288 = Icges(5,5) * t246 - Icges(5,6) * t245;
t242 = t269 * pkin(1) + t266 * pkin(7);
t287 = -V_base(4) * pkin(6) + t249 * t242 + V_base(2);
t164 = t269 * t309 + t207;
t286 = V_base(4) * t241 - t242 * V_base(5) + V_base(3);
t206 = V_base(5) + (-qJD(4) + t307) * t269;
t285 = pkin(11) * t245 + t246 * t332;
t284 = (-Icges(5,3) * t269 + t266 * t288) * t206 + (Icges(5,3) * t266 + t269 * t288) * t207 + (Icges(5,5) * t245 + Icges(5,6) * t246) * t249;
t283 = (-Icges(4,3) * t269 + t266 * t289) * t218 + (Icges(4,3) * t266 + t269 * t289) * t219 + (Icges(4,5) * t254 + Icges(4,6) * t256) * t249;
t282 = (-Icges(3,3) * t269 + t266 * t290) * t243 + (Icges(3,3) * t266 + t269 * t290) * t244 + (Icges(3,5) * t265 + Icges(3,6) * t268) * t249;
t163 = t266 * t309 + t206;
t177 = pkin(8) * t266 + t269 * t333;
t281 = t244 * t176 - t177 * t243 + t286;
t280 = t249 * t177 - t244 * t335 + t287;
t190 = t300 * t266;
t212 = pkin(4) * t245 - pkin(10) * t246;
t279 = t206 * t212 + (-t190 + t303) * t249 + t301;
t149 = pkin(9) * t266 + t269 * t311;
t278 = t219 * t148 - t149 * t218 + t281;
t277 = t249 * t149 - t219 * t334 + t280;
t191 = t300 * t269;
t276 = t207 * t190 - t191 * t206 + t278;
t275 = t249 * t191 - t207 * t212 + t277;
t170 = -Icges(5,6) * t269 + t266 * t291;
t171 = Icges(5,6) * t266 + t269 * t291;
t172 = -Icges(5,5) * t269 + t266 * t294;
t173 = Icges(5,5) * t266 + t269 * t294;
t209 = Icges(5,2) * t246 + t324;
t210 = Icges(5,1) * t245 + t323;
t274 = (-t171 * t245 + t173 * t246) * t207 + (-t170 * t245 + t172 * t246) * t206 + (-t209 * t245 + t210 * t246) * t249;
t180 = -Icges(4,6) * t269 + t266 * t292;
t181 = Icges(4,6) * t266 + t269 * t292;
t182 = -Icges(4,5) * t269 + t266 * t295;
t183 = Icges(4,5) * t266 + t269 * t295;
t215 = Icges(4,2) * t256 + t326;
t216 = Icges(4,1) * t254 + t325;
t273 = (-t181 * t254 + t183 * t256) * t219 + (-t180 * t254 + t182 * t256) * t218 + (-t215 * t254 + t216 * t256) * t249;
t195 = -Icges(3,6) * t269 + t266 * t293;
t196 = Icges(3,6) * t266 + t269 * t293;
t197 = -Icges(3,5) * t269 + t266 * t296;
t198 = Icges(3,5) * t266 + t269 * t296;
t232 = Icges(3,2) * t268 + t328;
t235 = Icges(3,1) * t265 + t327;
t272 = (-t196 * t265 + t198 * t268) * t244 + (-t195 * t265 + t197 * t268) * t243 + (-t232 * t265 + t235 * t268) * t249;
t257 = Icges(2,4) * t269;
t240 = rSges(2,1) * t269 - rSges(2,2) * t266;
t239 = rSges(2,1) * t266 + rSges(2,2) * t269;
t238 = rSges(3,1) * t265 + rSges(3,2) * t268;
t237 = Icges(2,1) * t269 - t329;
t236 = Icges(2,1) * t266 + t257;
t234 = -Icges(2,2) * t266 + t257;
t233 = Icges(2,2) * t269 + t329;
t226 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t225 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t224 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t220 = -qJD(5) * t246 + t249;
t217 = rSges(4,1) * t254 + rSges(4,2) * t256;
t211 = rSges(5,1) * t245 + rSges(5,2) * t246;
t205 = t246 * t313 + t316;
t204 = -t246 * t315 + t314;
t203 = t246 * t314 - t315;
t202 = -t246 * t316 - t313;
t200 = rSges(3,3) * t266 + t269 * t299;
t199 = -rSges(3,3) * t269 + t266 * t299;
t192 = (-qJD(5) - qJD(6)) * t246 + t249;
t189 = t246 * t317 + t320;
t188 = -t246 * t319 + t318;
t187 = t246 * t318 - t319;
t186 = -t246 * t320 - t317;
t185 = rSges(4,3) * t266 + t269 * t298;
t184 = -rSges(4,3) * t269 + t266 * t298;
t175 = rSges(5,3) * t266 + t269 * t297;
t174 = -rSges(5,3) * t269 + t266 * t297;
t166 = V_base(5) * rSges(2,3) - t239 * t249 + t306;
t165 = t240 * t249 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t162 = t239 * V_base(4) - t240 * V_base(5) + V_base(3);
t161 = -rSges(6,3) * t246 + (rSges(6,1) * t267 - rSges(6,2) * t264) * t245;
t160 = -Icges(6,5) * t246 + (Icges(6,1) * t267 - Icges(6,4) * t264) * t245;
t159 = -Icges(6,6) * t246 + (Icges(6,4) * t267 - Icges(6,2) * t264) * t245;
t158 = -Icges(6,3) * t246 + (Icges(6,5) * t267 - Icges(6,6) * t264) * t245;
t156 = -rSges(7,3) * t246 + (rSges(7,1) * t255 - rSges(7,2) * t253) * t245;
t154 = -Icges(7,5) * t246 + (Icges(7,1) * t255 - Icges(7,4) * t253) * t245;
t153 = -Icges(7,6) * t246 + (Icges(7,4) * t255 - Icges(7,2) * t253) * t245;
t152 = -Icges(7,3) * t246 + (Icges(7,5) * t255 - Icges(7,6) * t253) * t245;
t150 = -pkin(11) * t246 + t245 * t332;
t147 = t269 * t308 + t164;
t146 = t266 * t308 + t163;
t142 = rSges(6,1) * t205 + rSges(6,2) * t204 + rSges(6,3) * t321;
t141 = rSges(6,1) * t203 + rSges(6,2) * t202 + rSges(6,3) * t322;
t140 = Icges(6,1) * t205 + Icges(6,4) * t204 + Icges(6,5) * t321;
t139 = Icges(6,1) * t203 + Icges(6,4) * t202 + Icges(6,5) * t322;
t138 = Icges(6,4) * t205 + Icges(6,2) * t204 + Icges(6,6) * t321;
t137 = Icges(6,4) * t203 + Icges(6,2) * t202 + Icges(6,6) * t322;
t136 = Icges(6,5) * t205 + Icges(6,6) * t204 + Icges(6,3) * t321;
t135 = Icges(6,5) * t203 + Icges(6,6) * t202 + Icges(6,3) * t322;
t134 = pkin(5) * t316 + t269 * t285;
t133 = -pkin(5) * t315 + t266 * t285;
t132 = rSges(7,1) * t189 + rSges(7,2) * t188 + rSges(7,3) * t321;
t131 = rSges(7,1) * t187 + rSges(7,2) * t186 + rSges(7,3) * t322;
t130 = Icges(7,1) * t189 + Icges(7,4) * t188 + Icges(7,5) * t321;
t129 = Icges(7,1) * t187 + Icges(7,4) * t186 + Icges(7,5) * t322;
t128 = Icges(7,4) * t189 + Icges(7,2) * t188 + Icges(7,6) * t321;
t127 = Icges(7,4) * t187 + Icges(7,2) * t186 + Icges(7,6) * t322;
t126 = Icges(7,5) * t189 + Icges(7,6) * t188 + Icges(7,3) * t321;
t125 = Icges(7,5) * t187 + Icges(7,6) * t186 + Icges(7,3) * t322;
t124 = t238 * t243 + (-t199 - t241) * t249 + t306;
t123 = t200 * t249 - t238 * t244 + t287;
t122 = t199 * t244 - t200 * t243 + t286;
t121 = t217 * t218 + (-t184 + t312) * t249 + t302;
t120 = t185 * t249 - t217 * t219 + t280;
t119 = t184 * t219 - t185 * t218 + t281;
t118 = t206 * t211 + (-t174 + t303) * t249 + t301;
t117 = t175 * t249 - t207 * t211 + t277;
t116 = t174 * t207 - t175 * t206 + t278;
t115 = -t141 * t220 + t161 * t163 + t279;
t114 = t142 * t220 - t161 * t164 + t275;
t113 = t141 * t164 - t142 * t163 + t276;
t112 = -t131 * t192 - t133 * t220 + t146 * t156 + t150 * t163 + t279;
t111 = t132 * t192 + t134 * t220 - t147 * t156 - t150 * t164 + t275;
t110 = t131 * t147 - t132 * t146 + t133 * t164 - t134 * t163 + t276;
t1 = V_base(5) * t249 * (Icges(2,5) * t266 + Icges(2,6) * t269) + t163 * ((t136 * t322 + t138 * t202 + t140 * t203) * t164 + (t135 * t322 + t202 * t137 + t203 * t139) * t163 + (t158 * t322 + t159 * t202 + t160 * t203) * t220) / 0.2e1 + t146 * ((t126 * t322 + t128 * t186 + t130 * t187) * t147 + (t125 * t322 + t186 * t127 + t187 * t129) * t146 + (t152 * t322 + t153 * t186 + t154 * t187) * t192) / 0.2e1 + t164 * ((t136 * t321 + t204 * t138 + t205 * t140) * t164 + (t135 * t321 + t137 * t204 + t139 * t205) * t163 + (t158 * t321 + t159 * t204 + t160 * t205) * t220) / 0.2e1 + t147 * ((t126 * t321 + t188 * t128 + t189 * t130) * t147 + (t125 * t321 + t127 * t188 + t129 * t189) * t146 + (t152 * t321 + t153 * t188 + t154 * t189) * t192) / 0.2e1 + t244 * (t282 * t266 + t272 * t269) / 0.2e1 + t243 * (t272 * t266 - t282 * t269) / 0.2e1 + t219 * (t283 * t266 + t273 * t269) / 0.2e1 + t218 * (t273 * t266 - t283 * t269) / 0.2e1 + t207 * (t284 * t266 + t274 * t269) / 0.2e1 + t206 * (t274 * t266 - t284 * t269) / 0.2e1 + t220 * ((-t135 * t163 - t136 * t164 - t158 * t220) * t246 + ((-t138 * t264 + t140 * t267) * t164 + (-t137 * t264 + t139 * t267) * t163 + (-t159 * t264 + t160 * t267) * t220) * t245) / 0.2e1 + ((t196 * t268 + t198 * t265) * t244 + (t195 * t268 + t197 * t265) * t243 + (t181 * t256 + t183 * t254) * t219 + (t180 * t256 + t182 * t254) * t218 + (t171 * t246 + t173 * t245) * t207 + (t170 * t246 + t172 * t245) * t206 + (t246 * t209 + t245 * t210 + t256 * t215 + t254 * t216 + t268 * t232 + t265 * t235 + Icges(2,3)) * t249) * t249 / 0.2e1 + t192 * ((-t125 * t146 - t126 * t147 - t152 * t192) * t246 + ((-t128 * t253 + t130 * t255) * t147 + (-t127 * t253 + t129 * t255) * t146 + (-t153 * t253 + t154 * t255) * t192) * t245) / 0.2e1 + ((-t233 * t266 + t236 * t269 + Icges(1,4)) * V_base(5) + (-t266 * t234 + t269 * t237 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t269 * t233 + t266 * t236 + Icges(1,2)) * V_base(5) + (t234 * t269 + t237 * t266 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + m(1) * (t224 ^ 2 + t225 ^ 2 + t226 ^ 2) / 0.2e1 + m(2) * (t162 ^ 2 + t165 ^ 2 + t166 ^ 2) / 0.2e1 + m(3) * (t122 ^ 2 + t123 ^ 2 + t124 ^ 2) / 0.2e1 + m(6) * (t113 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + m(5) * (t116 ^ 2 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + m(4) * (t119 ^ 2 + t120 ^ 2 + t121 ^ 2) / 0.2e1 + m(7) * (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + t249 * V_base(4) * (Icges(2,5) * t269 - Icges(2,6) * t266);
T  = t1;
