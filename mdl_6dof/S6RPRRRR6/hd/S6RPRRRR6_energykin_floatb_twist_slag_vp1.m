% Calculate kinetic energy for
% S6RPRRRR6
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
% Datum: 2019-03-09 07:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRR6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR6_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR6_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRRR6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR6_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR6_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR6_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRR6_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:12:10
% EndTime: 2019-03-09 07:12:13
% DurationCPUTime: 3.59s
% Computational Cost: add. (2445->389), mult. (2352->586), div. (0->0), fcn. (2246->12), ass. (0->192)
t267 = sin(pkin(11));
t335 = pkin(2) * t267;
t268 = cos(pkin(11));
t333 = pkin(2) * t268;
t272 = cos(qJ(4));
t332 = t272 * pkin(4);
t271 = sin(qJ(1));
t330 = Icges(2,4) * t271;
t329 = Icges(3,4) * t267;
t328 = Icges(3,4) * t268;
t264 = pkin(11) + qJ(3);
t253 = sin(t264);
t327 = Icges(4,4) * t253;
t254 = cos(t264);
t326 = Icges(4,4) * t254;
t325 = t253 * t271;
t273 = cos(qJ(1));
t324 = t253 * t273;
t323 = t254 * t271;
t322 = t254 * t273;
t266 = qJ(4) + qJ(5);
t258 = sin(t266);
t321 = t258 * t271;
t320 = t258 * t273;
t259 = cos(t266);
t319 = t259 * t271;
t318 = t259 * t273;
t270 = sin(qJ(4));
t317 = t270 * t271;
t316 = t270 * t273;
t315 = t271 * t272;
t314 = t272 * t273;
t178 = -pkin(7) * t273 + t271 * t333;
t242 = pkin(1) * t271 - qJ(2) * t273;
t312 = -t178 - t242;
t311 = pkin(5) * t259;
t309 = qJD(4) * t253;
t308 = qJD(5) * t253;
t307 = qJD(6) * t253;
t306 = -qJD(4) - qJD(5);
t305 = V_base(4) * t242 + V_base(3);
t304 = V_base(5) * pkin(6) + V_base(1);
t247 = qJD(3) * t271 + V_base(4);
t255 = V_base(6) + qJD(1);
t301 = pkin(5) * t258;
t300 = qJD(2) * t271 + t304;
t210 = t273 * t309 + t247;
t299 = V_base(5) * t335 + t300;
t298 = pkin(3) * t254 + pkin(8) * t253;
t246 = -qJD(3) * t273 + V_base(5);
t297 = rSges(3,1) * t268 - rSges(3,2) * t267;
t296 = rSges(4,1) * t254 - rSges(4,2) * t253;
t170 = t273 * t308 + t210;
t295 = Icges(3,1) * t268 - t329;
t294 = Icges(4,1) * t254 - t327;
t293 = -Icges(3,2) * t267 + t328;
t292 = -Icges(4,2) * t253 + t326;
t291 = Icges(3,5) * t268 - Icges(3,6) * t267;
t290 = Icges(4,5) * t254 - Icges(4,6) * t253;
t244 = pkin(1) * t273 + qJ(2) * t271;
t289 = -qJD(2) * t273 + t255 * t244 + V_base(2);
t209 = t271 * t309 + t246;
t169 = t271 * t308 + t209;
t288 = pkin(9) * t253 + t254 * t332;
t287 = (-Icges(4,3) * t273 + t271 * t290) * t246 + (Icges(4,3) * t271 + t273 * t290) * t247 + (Icges(4,5) * t253 + Icges(4,6) * t254) * t255;
t286 = pkin(10) * t253 + t254 * t311;
t179 = pkin(7) * t271 + t273 * t333;
t285 = V_base(4) * t178 + (-t179 - t244) * V_base(5) + t305;
t284 = (-Icges(3,3) * t273 + t271 * t291) * V_base(5) + (Icges(3,3) * t271 + t273 * t291) * V_base(4) + (Icges(3,5) * t267 + Icges(3,6) * t268) * t255;
t206 = t298 * t271;
t220 = t253 * pkin(3) - t254 * pkin(8);
t283 = t246 * t220 + (-t206 + t312) * t255 + t299;
t207 = t298 * t273;
t282 = t247 * t206 - t207 * t246 + t285;
t281 = t255 * t179 + (-pkin(6) - t335) * V_base(4) + t289;
t143 = -pkin(4) * t316 + t271 * t288;
t156 = -pkin(9) * t254 + t253 * t332;
t223 = -qJD(4) * t254 + t255;
t280 = -t143 * t223 + t209 * t156 + t283;
t144 = pkin(4) * t317 + t273 * t288;
t279 = t210 * t143 - t144 * t209 + t282;
t278 = t255 * t207 - t220 * t247 + t281;
t277 = t223 * t144 - t156 * t210 + t278;
t183 = -Icges(4,6) * t273 + t271 * t292;
t184 = Icges(4,6) * t271 + t273 * t292;
t185 = -Icges(4,5) * t273 + t271 * t294;
t186 = Icges(4,5) * t271 + t273 * t294;
t217 = Icges(4,2) * t254 + t327;
t218 = Icges(4,1) * t253 + t326;
t276 = (-t184 * t253 + t186 * t254) * t247 + (-t183 * t253 + t185 * t254) * t246 + (-t217 * t253 + t218 * t254) * t255;
t200 = -Icges(3,6) * t273 + t271 * t293;
t201 = Icges(3,6) * t271 + t273 * t293;
t202 = -Icges(3,5) * t273 + t271 * t295;
t203 = Icges(3,5) * t271 + t273 * t295;
t229 = Icges(3,2) * t268 + t329;
t230 = Icges(3,1) * t267 + t328;
t275 = (-t201 * t267 + t203 * t268) * V_base(4) + (-t200 * t267 + t202 * t268) * V_base(5) + (-t229 * t267 + t230 * t268) * t255;
t261 = qJ(6) + t266;
t260 = Icges(2,4) * t273;
t251 = cos(t261);
t250 = sin(t261);
t245 = rSges(2,1) * t273 - rSges(2,2) * t271;
t243 = rSges(2,1) * t271 + rSges(2,2) * t273;
t237 = Icges(2,1) * t273 - t330;
t236 = Icges(2,1) * t271 + t260;
t235 = -Icges(2,2) * t271 + t260;
t234 = Icges(2,2) * t273 + t330;
t233 = Icges(2,5) * t273 - Icges(2,6) * t271;
t232 = Icges(2,5) * t271 + Icges(2,6) * t273;
t231 = rSges(3,1) * t267 + rSges(3,2) * t268;
t226 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t225 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t224 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t219 = rSges(4,1) * t253 + rSges(4,2) * t254;
t214 = t254 * t314 + t317;
t213 = -t254 * t316 + t315;
t212 = t254 * t315 - t316;
t211 = -t254 * t317 - t314;
t208 = t254 * t306 + t255;
t205 = rSges(3,3) * t271 + t273 * t297;
t204 = -rSges(3,3) * t273 + t271 * t297;
t197 = t254 * t318 + t321;
t196 = -t254 * t320 + t319;
t195 = t254 * t319 - t320;
t194 = -t254 * t321 - t318;
t192 = rSges(4,3) * t271 + t273 * t296;
t191 = -rSges(4,3) * t273 + t271 * t296;
t190 = t250 * t271 + t251 * t322;
t189 = -t250 * t322 + t251 * t271;
t188 = -t250 * t273 + t251 * t323;
t187 = -t250 * t323 - t251 * t273;
t177 = (-qJD(6) + t306) * t254 + t255;
t176 = -rSges(5,3) * t254 + (rSges(5,1) * t272 - rSges(5,2) * t270) * t253;
t175 = V_base(5) * rSges(2,3) - t243 * t255 + t304;
t174 = t245 * t255 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t173 = -Icges(5,5) * t254 + (Icges(5,1) * t272 - Icges(5,4) * t270) * t253;
t172 = -Icges(5,6) * t254 + (Icges(5,4) * t272 - Icges(5,2) * t270) * t253;
t171 = -Icges(5,3) * t254 + (Icges(5,5) * t272 - Icges(5,6) * t270) * t253;
t167 = t243 * V_base(4) - t245 * V_base(5) + V_base(3);
t165 = -rSges(6,3) * t254 + (rSges(6,1) * t259 - rSges(6,2) * t258) * t253;
t164 = -Icges(6,5) * t254 + (Icges(6,1) * t259 - Icges(6,4) * t258) * t253;
t163 = -Icges(6,6) * t254 + (Icges(6,4) * t259 - Icges(6,2) * t258) * t253;
t162 = -Icges(6,3) * t254 + (Icges(6,5) * t259 - Icges(6,6) * t258) * t253;
t160 = -rSges(7,3) * t254 + (rSges(7,1) * t251 - rSges(7,2) * t250) * t253;
t159 = -Icges(7,5) * t254 + (Icges(7,1) * t251 - Icges(7,4) * t250) * t253;
t158 = -Icges(7,6) * t254 + (Icges(7,4) * t251 - Icges(7,2) * t250) * t253;
t157 = -Icges(7,3) * t254 + (Icges(7,5) * t251 - Icges(7,6) * t250) * t253;
t155 = t273 * t307 + t170;
t154 = t271 * t307 + t169;
t153 = -pkin(10) * t254 + t253 * t311;
t152 = rSges(5,1) * t214 + rSges(5,2) * t213 + rSges(5,3) * t324;
t151 = rSges(5,1) * t212 + rSges(5,2) * t211 + rSges(5,3) * t325;
t150 = Icges(5,1) * t214 + Icges(5,4) * t213 + Icges(5,5) * t324;
t149 = Icges(5,1) * t212 + Icges(5,4) * t211 + Icges(5,5) * t325;
t148 = Icges(5,4) * t214 + Icges(5,2) * t213 + Icges(5,6) * t324;
t147 = Icges(5,4) * t212 + Icges(5,2) * t211 + Icges(5,6) * t325;
t146 = Icges(5,5) * t214 + Icges(5,6) * t213 + Icges(5,3) * t324;
t145 = Icges(5,5) * t212 + Icges(5,6) * t211 + Icges(5,3) * t325;
t141 = rSges(6,1) * t197 + rSges(6,2) * t196 + rSges(6,3) * t324;
t140 = rSges(6,1) * t195 + rSges(6,2) * t194 + rSges(6,3) * t325;
t139 = Icges(6,1) * t197 + Icges(6,4) * t196 + Icges(6,5) * t324;
t138 = Icges(6,1) * t195 + Icges(6,4) * t194 + Icges(6,5) * t325;
t137 = Icges(6,4) * t197 + Icges(6,2) * t196 + Icges(6,6) * t324;
t136 = Icges(6,4) * t195 + Icges(6,2) * t194 + Icges(6,6) * t325;
t135 = Icges(6,5) * t197 + Icges(6,6) * t196 + Icges(6,3) * t324;
t134 = Icges(6,5) * t195 + Icges(6,6) * t194 + Icges(6,3) * t325;
t133 = rSges(7,1) * t190 + rSges(7,2) * t189 + rSges(7,3) * t324;
t132 = rSges(7,1) * t188 + rSges(7,2) * t187 + rSges(7,3) * t325;
t131 = Icges(7,1) * t190 + Icges(7,4) * t189 + Icges(7,5) * t324;
t130 = Icges(7,1) * t188 + Icges(7,4) * t187 + Icges(7,5) * t325;
t129 = Icges(7,4) * t190 + Icges(7,2) * t189 + Icges(7,6) * t324;
t128 = Icges(7,4) * t188 + Icges(7,2) * t187 + Icges(7,6) * t325;
t127 = Icges(7,5) * t190 + Icges(7,6) * t189 + Icges(7,3) * t324;
t126 = Icges(7,5) * t188 + Icges(7,6) * t187 + Icges(7,3) * t325;
t124 = t231 * V_base(5) + (-t204 - t242) * t255 + t300;
t123 = t205 * t255 + (-pkin(6) - t231) * V_base(4) + t289;
t121 = t204 * V_base(4) + (-t205 - t244) * V_base(5) + t305;
t120 = t271 * t301 + t273 * t286;
t119 = t271 * t286 - t273 * t301;
t118 = t219 * t246 + (-t191 + t312) * t255 + t299;
t117 = t192 * t255 - t219 * t247 + t281;
t116 = t191 * t247 - t192 * t246 + t285;
t115 = -t151 * t223 + t176 * t209 + t283;
t114 = t152 * t223 - t176 * t210 + t278;
t113 = t151 * t210 - t152 * t209 + t282;
t112 = -t140 * t208 + t165 * t169 + t280;
t111 = t141 * t208 - t165 * t170 + t277;
t110 = t140 * t170 - t141 * t169 + t279;
t109 = -t119 * t208 - t132 * t177 + t153 * t169 + t154 * t160 + t280;
t108 = t120 * t208 + t133 * t177 - t153 * t170 - t155 * t160 + t277;
t107 = t119 * t170 - t120 * t169 + t132 * t155 - t133 * t154 + t279;
t1 = m(3) * (t121 ^ 2 + t123 ^ 2 + t124 ^ 2) / 0.2e1 + m(5) * (t113 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + m(4) * (t116 ^ 2 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + m(7) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + m(6) * (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + t247 * (t287 * t271 + t276 * t273) / 0.2e1 + t246 * (t276 * t271 - t287 * t273) / 0.2e1 + m(1) * (t224 ^ 2 + t225 ^ 2 + t226 ^ 2) / 0.2e1 + m(2) * (t167 ^ 2 + t174 ^ 2 + t175 ^ 2) / 0.2e1 + t210 * ((t146 * t324 + t213 * t148 + t214 * t150) * t210 + (t145 * t324 + t147 * t213 + t149 * t214) * t209 + (t171 * t324 + t172 * t213 + t173 * t214) * t223) / 0.2e1 + t170 * ((t135 * t324 + t196 * t137 + t197 * t139) * t170 + (t134 * t324 + t136 * t196 + t138 * t197) * t169 + (t162 * t324 + t163 * t196 + t164 * t197) * t208) / 0.2e1 + t155 * ((t127 * t324 + t189 * t129 + t190 * t131) * t155 + (t126 * t324 + t128 * t189 + t130 * t190) * t154 + (t157 * t324 + t158 * t189 + t159 * t190) * t177) / 0.2e1 + t209 * ((t146 * t325 + t148 * t211 + t150 * t212) * t210 + (t145 * t325 + t211 * t147 + t212 * t149) * t209 + (t171 * t325 + t172 * t211 + t173 * t212) * t223) / 0.2e1 + t169 * ((t135 * t325 + t137 * t194 + t139 * t195) * t170 + (t134 * t325 + t194 * t136 + t195 * t138) * t169 + (t162 * t325 + t163 * t194 + t164 * t195) * t208) / 0.2e1 + t154 * ((t127 * t325 + t129 * t187 + t131 * t188) * t155 + (t126 * t325 + t187 * t128 + t188 * t130) * t154 + (t157 * t325 + t158 * t187 + t159 * t188) * t177) / 0.2e1 + t177 * ((-t126 * t154 - t127 * t155 - t157 * t177) * t254 + ((-t129 * t250 + t131 * t251) * t155 + (-t128 * t250 + t130 * t251) * t154 + (-t158 * t250 + t159 * t251) * t177) * t253) / 0.2e1 + t208 * ((-t134 * t169 - t135 * t170 - t162 * t208) * t254 + ((-t137 * t258 + t139 * t259) * t170 + (-t136 * t258 + t138 * t259) * t169 + (-t163 * t258 + t164 * t259) * t208) * t253) / 0.2e1 + t223 * ((-t145 * t209 - t146 * t210 - t171 * t223) * t254 + ((-t148 * t270 + t150 * t272) * t210 + (-t147 * t270 + t149 * t272) * t209 + (-t172 * t270 + t173 * t272) * t223) * t253) / 0.2e1 + ((t184 * t254 + t186 * t253) * t247 + (t183 * t254 + t185 * t253) * t246 + (t200 * t268 + t202 * t267 + t232) * V_base(5) + (t201 * t268 + t203 * t267 + t233) * V_base(4) + (t217 * t254 + t218 * t253 + t229 * t268 + t230 * t267 + Icges(2,3)) * t255) * t255 / 0.2e1 + (t233 * t255 + t284 * t271 + t275 * t273 + (-t234 * t271 + t236 * t273 + Icges(1,4)) * V_base(5) + (-t235 * t271 + t237 * t273 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t232 * t255 + t275 * t271 - t284 * t273 + (t234 * t273 + t236 * t271 + Icges(1,2)) * V_base(5) + (t235 * t273 + t237 * t271 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
