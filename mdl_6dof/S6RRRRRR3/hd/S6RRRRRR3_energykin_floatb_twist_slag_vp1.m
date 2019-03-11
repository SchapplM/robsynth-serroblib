% Calculate kinetic energy for
% S6RRRRRR3
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
% Datum: 2019-03-10 03:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR3_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR3_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRRR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR3_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR3_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR3_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRR3_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:39:07
% EndTime: 2019-03-10 03:39:10
% DurationCPUTime: 3.42s
% Computational Cost: add. (2535->386), mult. (2442->594), div. (0->0), fcn. (2336->12), ass. (0->187)
t266 = sin(qJ(2));
t328 = pkin(2) * t266;
t269 = cos(qJ(2));
t326 = pkin(2) * t269;
t268 = cos(qJ(4));
t325 = t268 * pkin(4);
t267 = sin(qJ(1));
t322 = Icges(2,4) * t267;
t321 = Icges(3,4) * t266;
t320 = Icges(3,4) * t269;
t264 = qJ(2) + qJ(3);
t255 = sin(t264);
t319 = Icges(4,4) * t255;
t257 = cos(t264);
t318 = Icges(4,4) * t257;
t317 = t255 * t267;
t270 = cos(qJ(1));
t316 = t255 * t270;
t315 = t257 * t267;
t314 = t257 * t270;
t265 = sin(qJ(4));
t313 = t265 * t267;
t312 = t265 * t270;
t311 = t267 * t268;
t310 = t268 * t270;
t263 = qJ(4) + qJ(5);
t175 = -pkin(8) * t270 + t267 * t326;
t243 = t267 * pkin(1) - t270 * pkin(7);
t309 = -t175 - t243;
t256 = cos(t263);
t308 = pkin(5) * t256;
t306 = qJD(4) * t255;
t305 = qJD(5) * t255;
t304 = qJD(6) * t255;
t303 = -qJD(4) - qJD(5);
t302 = V_base(5) * pkin(6) + V_base(1);
t246 = qJD(2) * t267 + V_base(4);
t251 = V_base(6) + qJD(1);
t254 = sin(t263);
t299 = pkin(5) * t254;
t245 = -qJD(2) * t270 + V_base(5);
t298 = t245 * t328 + t302;
t218 = qJD(3) * t267 + t246;
t297 = pkin(3) * t257 + pkin(9) * t255;
t296 = rSges(3,1) * t269 - rSges(3,2) * t266;
t295 = rSges(4,1) * t257 - rSges(4,2) * t255;
t191 = t270 * t306 + t218;
t294 = Icges(3,1) * t269 - t321;
t293 = Icges(4,1) * t257 - t319;
t292 = -Icges(3,2) * t266 + t320;
t291 = -Icges(4,2) * t255 + t318;
t290 = Icges(3,5) * t269 - Icges(3,6) * t266;
t289 = Icges(4,5) * t257 - Icges(4,6) * t255;
t244 = t270 * pkin(1) + t267 * pkin(7);
t288 = -V_base(4) * pkin(6) + t251 * t244 + V_base(2);
t159 = t270 * t305 + t191;
t287 = V_base(4) * t243 - t244 * V_base(5) + V_base(3);
t217 = V_base(5) + (-qJD(2) - qJD(3)) * t270;
t190 = t267 * t306 + t217;
t286 = pkin(10) * t255 + t257 * t325;
t285 = (-Icges(4,3) * t270 + t267 * t289) * t217 + (Icges(4,3) * t267 + t270 * t289) * t218 + (Icges(4,5) * t255 + Icges(4,6) * t257) * t251;
t284 = (-Icges(3,3) * t270 + t267 * t290) * t245 + (Icges(3,3) * t267 + t270 * t290) * t246 + (Icges(3,5) * t266 + Icges(3,6) * t269) * t251;
t283 = pkin(11) * t255 + t257 * t308;
t158 = t267 * t305 + t190;
t204 = t297 * t267;
t216 = pkin(3) * t255 - pkin(9) * t257;
t282 = t217 * t216 + (-t204 + t309) * t251 + t298;
t176 = pkin(8) * t267 + t270 * t326;
t281 = t246 * t175 - t176 * t245 + t287;
t280 = t251 * t176 - t246 * t328 + t288;
t139 = -pkin(4) * t312 + t267 * t286;
t153 = -pkin(10) * t257 + t255 * t325;
t225 = -qJD(4) * t257 + t251;
t279 = -t139 * t225 + t190 * t153 + t282;
t205 = t297 * t270;
t278 = t218 * t204 - t205 * t217 + t281;
t277 = t251 * t205 - t216 * t218 + t280;
t140 = pkin(4) * t313 + t270 * t286;
t276 = t191 * t139 - t140 * t190 + t278;
t275 = t225 * t140 - t153 * t191 + t277;
t180 = -Icges(4,6) * t270 + t267 * t291;
t181 = Icges(4,6) * t267 + t270 * t291;
t182 = -Icges(4,5) * t270 + t267 * t293;
t183 = Icges(4,5) * t267 + t270 * t293;
t213 = Icges(4,2) * t257 + t319;
t214 = Icges(4,1) * t255 + t318;
t274 = (-t181 * t255 + t183 * t257) * t218 + (-t180 * t255 + t182 * t257) * t217 + (-t213 * t255 + t214 * t257) * t251;
t198 = -Icges(3,6) * t270 + t267 * t292;
t199 = Icges(3,6) * t267 + t270 * t292;
t200 = -Icges(3,5) * t270 + t267 * t294;
t201 = Icges(3,5) * t267 + t270 * t294;
t230 = Icges(3,2) * t269 + t321;
t233 = Icges(3,1) * t266 + t320;
t273 = (-t199 * t266 + t201 * t269) * t246 + (-t198 * t266 + t200 * t269) * t245 + (-t230 * t266 + t233 * t269) * t251;
t259 = qJ(6) + t263;
t258 = Icges(2,4) * t270;
t248 = cos(t259);
t247 = sin(t259);
t238 = rSges(2,1) * t270 - rSges(2,2) * t267;
t237 = rSges(2,1) * t267 + rSges(2,2) * t270;
t236 = rSges(3,1) * t266 + rSges(3,2) * t269;
t235 = Icges(2,1) * t270 - t322;
t234 = Icges(2,1) * t267 + t258;
t232 = -Icges(2,2) * t267 + t258;
t231 = Icges(2,2) * t270 + t322;
t224 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t223 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t222 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t215 = rSges(4,1) * t255 + rSges(4,2) * t257;
t210 = t257 * t310 + t313;
t209 = -t257 * t312 + t311;
t208 = t257 * t311 - t312;
t207 = -t257 * t313 - t310;
t206 = t257 * t303 + t251;
t203 = rSges(3,3) * t267 + t270 * t296;
t202 = -rSges(3,3) * t270 + t267 * t296;
t195 = t254 * t267 + t256 * t314;
t194 = -t254 * t314 + t256 * t267;
t193 = -t254 * t270 + t256 * t315;
t192 = -t254 * t315 - t256 * t270;
t189 = rSges(4,3) * t267 + t270 * t295;
t188 = -rSges(4,3) * t270 + t267 * t295;
t187 = t247 * t267 + t248 * t314;
t186 = -t247 * t314 + t248 * t267;
t185 = -t247 * t270 + t248 * t315;
t184 = -t247 * t315 - t248 * t270;
t174 = (-qJD(6) + t303) * t257 + t251;
t173 = -rSges(5,3) * t257 + (rSges(5,1) * t268 - rSges(5,2) * t265) * t255;
t172 = -Icges(5,5) * t257 + (Icges(5,1) * t268 - Icges(5,4) * t265) * t255;
t171 = -Icges(5,6) * t257 + (Icges(5,4) * t268 - Icges(5,2) * t265) * t255;
t170 = -Icges(5,3) * t257 + (Icges(5,5) * t268 - Icges(5,6) * t265) * t255;
t169 = V_base(5) * rSges(2,3) - t237 * t251 + t302;
t168 = t238 * t251 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t166 = t237 * V_base(4) - t238 * V_base(5) + V_base(3);
t165 = -rSges(6,3) * t257 + (rSges(6,1) * t256 - rSges(6,2) * t254) * t255;
t164 = -Icges(6,5) * t257 + (Icges(6,1) * t256 - Icges(6,4) * t254) * t255;
t163 = -Icges(6,6) * t257 + (Icges(6,4) * t256 - Icges(6,2) * t254) * t255;
t162 = -Icges(6,3) * t257 + (Icges(6,5) * t256 - Icges(6,6) * t254) * t255;
t160 = -rSges(7,3) * t257 + (rSges(7,1) * t248 - rSges(7,2) * t247) * t255;
t156 = -Icges(7,5) * t257 + (Icges(7,1) * t248 - Icges(7,4) * t247) * t255;
t155 = -Icges(7,6) * t257 + (Icges(7,4) * t248 - Icges(7,2) * t247) * t255;
t154 = -Icges(7,3) * t257 + (Icges(7,5) * t248 - Icges(7,6) * t247) * t255;
t151 = t270 * t304 + t159;
t150 = t267 * t304 + t158;
t149 = -pkin(11) * t257 + t255 * t308;
t148 = rSges(5,1) * t210 + rSges(5,2) * t209 + rSges(5,3) * t316;
t147 = rSges(5,1) * t208 + rSges(5,2) * t207 + rSges(5,3) * t317;
t146 = Icges(5,1) * t210 + Icges(5,4) * t209 + Icges(5,5) * t316;
t145 = Icges(5,1) * t208 + Icges(5,4) * t207 + Icges(5,5) * t317;
t144 = Icges(5,4) * t210 + Icges(5,2) * t209 + Icges(5,6) * t316;
t143 = Icges(5,4) * t208 + Icges(5,2) * t207 + Icges(5,6) * t317;
t142 = Icges(5,5) * t210 + Icges(5,6) * t209 + Icges(5,3) * t316;
t141 = Icges(5,5) * t208 + Icges(5,6) * t207 + Icges(5,3) * t317;
t137 = rSges(6,1) * t195 + rSges(6,2) * t194 + rSges(6,3) * t316;
t136 = rSges(6,1) * t193 + rSges(6,2) * t192 + rSges(6,3) * t317;
t135 = Icges(6,1) * t195 + Icges(6,4) * t194 + Icges(6,5) * t316;
t134 = Icges(6,1) * t193 + Icges(6,4) * t192 + Icges(6,5) * t317;
t133 = Icges(6,4) * t195 + Icges(6,2) * t194 + Icges(6,6) * t316;
t132 = Icges(6,4) * t193 + Icges(6,2) * t192 + Icges(6,6) * t317;
t131 = Icges(6,5) * t195 + Icges(6,6) * t194 + Icges(6,3) * t316;
t130 = Icges(6,5) * t193 + Icges(6,6) * t192 + Icges(6,3) * t317;
t129 = rSges(7,1) * t187 + rSges(7,2) * t186 + rSges(7,3) * t316;
t128 = rSges(7,1) * t185 + rSges(7,2) * t184 + rSges(7,3) * t317;
t127 = Icges(7,1) * t187 + Icges(7,4) * t186 + Icges(7,5) * t316;
t126 = Icges(7,1) * t185 + Icges(7,4) * t184 + Icges(7,5) * t317;
t125 = Icges(7,4) * t187 + Icges(7,2) * t186 + Icges(7,6) * t316;
t124 = Icges(7,4) * t185 + Icges(7,2) * t184 + Icges(7,6) * t317;
t123 = Icges(7,5) * t187 + Icges(7,6) * t186 + Icges(7,3) * t316;
t122 = Icges(7,5) * t185 + Icges(7,6) * t184 + Icges(7,3) * t317;
t120 = t236 * t245 + (-t202 - t243) * t251 + t302;
t119 = t203 * t251 - t236 * t246 + t288;
t117 = t267 * t299 + t270 * t283;
t116 = t267 * t283 - t270 * t299;
t115 = t202 * t246 - t203 * t245 + t287;
t114 = t215 * t217 + (-t188 + t309) * t251 + t298;
t113 = t189 * t251 - t215 * t218 + t280;
t112 = t188 * t218 - t189 * t217 + t281;
t111 = -t147 * t225 + t173 * t190 + t282;
t110 = t148 * t225 - t173 * t191 + t277;
t109 = t147 * t191 - t148 * t190 + t278;
t108 = -t136 * t206 + t158 * t165 + t279;
t107 = t137 * t206 - t159 * t165 + t275;
t106 = t136 * t159 - t137 * t158 + t276;
t105 = -t116 * t206 - t128 * t174 + t149 * t158 + t150 * t160 + t279;
t104 = t117 * t206 + t129 * t174 - t149 * t159 - t151 * t160 + t275;
t103 = t116 * t159 - t117 * t158 + t128 * t151 - t129 * t150 + t276;
t1 = t174 * ((-t122 * t150 - t123 * t151 - t154 * t174) * t257 + ((-t125 * t247 + t127 * t248) * t151 + (-t124 * t247 + t126 * t248) * t150 + (-t155 * t247 + t156 * t248) * t174) * t255) / 0.2e1 + m(1) * (t222 ^ 2 + t223 ^ 2 + t224 ^ 2) / 0.2e1 + m(2) * (t166 ^ 2 + t168 ^ 2 + t169 ^ 2) / 0.2e1 + m(3) * (t115 ^ 2 + t119 ^ 2 + t120 ^ 2) / 0.2e1 + m(5) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + m(4) * (t112 ^ 2 + t113 ^ 2 + t114 ^ 2) / 0.2e1 + m(7) * (t103 ^ 2 + t104 ^ 2 + t105 ^ 2) / 0.2e1 + m(6) * (t106 ^ 2 + t107 ^ 2 + t108 ^ 2) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + V_base(5) * t251 * (Icges(2,5) * t267 + Icges(2,6) * t270) + ((-t231 * t267 + t234 * t270 + Icges(1,4)) * V_base(5) + (-t267 * t232 + t270 * t235 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + t150 * ((t123 * t317 + t125 * t184 + t127 * t185) * t151 + (t122 * t317 + t184 * t124 + t126 * t185) * t150 + (t154 * t317 + t155 * t184 + t156 * t185) * t174) / 0.2e1 + t191 * ((t142 * t316 + t209 * t144 + t210 * t146) * t191 + (t141 * t316 + t143 * t209 + t145 * t210) * t190 + (t170 * t316 + t171 * t209 + t172 * t210) * t225) / 0.2e1 + t159 * ((t131 * t316 + t133 * t194 + t195 * t135) * t159 + (t130 * t316 + t132 * t194 + t134 * t195) * t158 + (t162 * t316 + t163 * t194 + t164 * t195) * t206) / 0.2e1 + t151 * ((t123 * t316 + t125 * t186 + t187 * t127) * t151 + (t122 * t316 + t124 * t186 + t126 * t187) * t150 + (t154 * t316 + t155 * t186 + t156 * t187) * t174) / 0.2e1 + t190 * ((t142 * t317 + t144 * t207 + t146 * t208) * t191 + (t141 * t317 + t207 * t143 + t145 * t208) * t190 + (t170 * t317 + t171 * t207 + t172 * t208) * t225) / 0.2e1 + t158 * ((t131 * t317 + t133 * t192 + t135 * t193) * t159 + (t130 * t317 + t192 * t132 + t193 * t134) * t158 + (t162 * t317 + t163 * t192 + t164 * t193) * t206) / 0.2e1 + t225 * ((-t141 * t190 - t142 * t191 - t170 * t225) * t257 + ((-t144 * t265 + t146 * t268) * t191 + (-t143 * t265 + t145 * t268) * t190 + (-t171 * t265 + t172 * t268) * t225) * t255) / 0.2e1 + ((t270 * t231 + t267 * t234 + Icges(1,2)) * V_base(5) + (t232 * t270 + t235 * t267 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + t251 * V_base(4) * (Icges(2,5) * t270 - Icges(2,6) * t267) + t246 * (t267 * t284 + t270 * t273) / 0.2e1 + t245 * (t267 * t273 - t270 * t284) / 0.2e1 + t218 * (t267 * t285 + t270 * t274) / 0.2e1 + t217 * (t267 * t274 - t270 * t285) / 0.2e1 + t206 * ((-t130 * t158 - t131 * t159 - t162 * t206) * t257 + ((-t133 * t254 + t135 * t256) * t159 + (-t132 * t254 + t134 * t256) * t158 + (-t163 * t254 + t164 * t256) * t206) * t255) / 0.2e1 + ((t199 * t269 + t201 * t266) * t246 + (t198 * t269 + t200 * t266) * t245 + (t181 * t257 + t183 * t255) * t218 + (t180 * t257 + t182 * t255) * t217 + (t257 * t213 + t255 * t214 + t269 * t230 + t266 * t233 + Icges(2,3)) * t251) * t251 / 0.2e1;
T  = t1;
