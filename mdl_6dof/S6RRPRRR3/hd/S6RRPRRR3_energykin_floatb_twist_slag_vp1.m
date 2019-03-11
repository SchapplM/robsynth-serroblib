% Calculate kinetic energy for
% S6RRPRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 13:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR3_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR3_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRRR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR3_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR3_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR3_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR3_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:21:40
% EndTime: 2019-03-09 13:21:44
% DurationCPUTime: 4.10s
% Computational Cost: add. (2471->383), mult. (2378->572), div. (0->0), fcn. (2272->12), ass. (0->187)
t341 = Icges(3,3) + Icges(4,3);
t264 = qJ(2) + pkin(11);
t253 = sin(t264);
t254 = cos(t264);
t269 = sin(qJ(2));
t272 = cos(qJ(2));
t340 = Icges(3,5) * t272 + Icges(4,5) * t254 - Icges(3,6) * t269 - Icges(4,6) * t253;
t270 = sin(qJ(1));
t273 = cos(qJ(1));
t324 = Icges(4,4) * t254;
t293 = -Icges(4,2) * t253 + t324;
t183 = -Icges(4,6) * t273 + t270 * t293;
t184 = Icges(4,6) * t270 + t273 * t293;
t325 = Icges(4,4) * t253;
t295 = Icges(4,1) * t254 - t325;
t185 = -Icges(4,5) * t273 + t270 * t295;
t186 = Icges(4,5) * t270 + t273 * t295;
t326 = Icges(3,4) * t272;
t294 = -Icges(3,2) * t269 + t326;
t200 = -Icges(3,6) * t273 + t270 * t294;
t201 = Icges(3,6) * t270 + t273 * t294;
t327 = Icges(3,4) * t269;
t296 = Icges(3,1) * t272 - t327;
t202 = -Icges(3,5) * t273 + t270 * t296;
t203 = Icges(3,5) * t270 + t273 * t296;
t217 = Icges(4,2) * t254 + t325;
t218 = Icges(4,1) * t253 + t324;
t232 = Icges(3,2) * t272 + t327;
t235 = Icges(3,1) * t269 + t326;
t247 = -qJD(2) * t273 + V_base(5);
t248 = qJD(2) * t270 + V_base(4);
t255 = V_base(6) + qJD(1);
t339 = (-t217 * t253 + t218 * t254 - t232 * t269 + t235 * t272) * t255 + (-t184 * t253 + t186 * t254 - t201 * t269 + t203 * t272) * t248 + (-t183 * t253 + t185 * t254 - t200 * t269 + t202 * t272) * t247;
t338 = (Icges(3,5) * t269 + Icges(4,5) * t253 + Icges(3,6) * t272 + Icges(4,6) * t254) * t255 + (t270 * t341 + t340 * t273) * t248 + (t340 * t270 - t273 * t341) * t247;
t334 = pkin(2) * t269;
t332 = pkin(2) * t272;
t271 = cos(qJ(4));
t331 = t271 * pkin(4);
t328 = Icges(2,4) * t270;
t323 = t253 * t270;
t322 = t253 * t273;
t321 = t254 * t270;
t320 = t254 * t273;
t266 = qJ(4) + qJ(5);
t258 = sin(t266);
t319 = t258 * t270;
t318 = t258 * t273;
t259 = cos(t266);
t317 = t259 * t270;
t316 = t259 * t273;
t268 = sin(qJ(4));
t315 = t268 * t270;
t314 = t268 * t273;
t313 = t270 * t271;
t312 = t271 * t273;
t178 = -qJ(3) * t273 + t270 * t332;
t245 = pkin(1) * t270 - pkin(7) * t273;
t311 = -t178 - t245;
t310 = pkin(5) * t259;
t308 = qJD(4) * t253;
t307 = qJD(5) * t253;
t306 = qJD(6) * t253;
t305 = -qJD(4) - qJD(5);
t304 = V_base(5) * pkin(6) + V_base(1);
t301 = pkin(5) * t258;
t210 = t273 * t308 + t248;
t300 = qJD(3) * t270 + t247 * t334 + t304;
t299 = pkin(3) * t254 + pkin(8) * t253;
t298 = rSges(3,1) * t272 - rSges(3,2) * t269;
t297 = rSges(4,1) * t254 - rSges(4,2) * t253;
t170 = t273 * t307 + t210;
t209 = t270 * t308 + t247;
t246 = pkin(1) * t273 + pkin(7) * t270;
t290 = -V_base(4) * pkin(6) + t255 * t246 + V_base(2);
t289 = V_base(4) * t245 - t246 * V_base(5) + V_base(3);
t169 = t270 * t307 + t209;
t288 = t248 * t178 + t289;
t287 = pkin(9) * t253 + t254 * t331;
t284 = pkin(10) * t253 + t254 * t310;
t179 = qJ(3) * t270 + t273 * t332;
t283 = -qJD(3) * t273 + t255 * t179 + t290;
t206 = t299 * t270;
t220 = t253 * pkin(3) - t254 * pkin(8);
t282 = t247 * t220 + (-t206 + t311) * t255 + t300;
t207 = t299 * t273;
t281 = t248 * t206 + (-t179 - t207) * t247 + t288;
t143 = -pkin(4) * t314 + t270 * t287;
t156 = -pkin(9) * t254 + t253 * t331;
t224 = -qJD(4) * t254 + t255;
t280 = -t143 * t224 + t209 * t156 + t282;
t144 = pkin(4) * t315 + t273 * t287;
t279 = t210 * t143 - t144 * t209 + t281;
t278 = t255 * t207 + (-t220 - t334) * t248 + t283;
t277 = t224 * t144 - t156 * t210 + t278;
t261 = qJ(6) + t266;
t260 = Icges(2,4) * t273;
t250 = cos(t261);
t249 = sin(t261);
t244 = rSges(2,1) * t273 - rSges(2,2) * t270;
t243 = rSges(2,1) * t270 + rSges(2,2) * t273;
t242 = rSges(3,1) * t269 + rSges(3,2) * t272;
t237 = Icges(2,1) * t273 - t328;
t236 = Icges(2,1) * t270 + t260;
t234 = -Icges(2,2) * t270 + t260;
t233 = Icges(2,2) * t273 + t328;
t227 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t226 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t225 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t219 = rSges(4,1) * t253 + rSges(4,2) * t254;
t214 = t254 * t312 + t315;
t213 = -t254 * t314 + t313;
t212 = t254 * t313 - t314;
t211 = -t254 * t315 - t312;
t208 = t254 * t305 + t255;
t205 = rSges(3,3) * t270 + t273 * t298;
t204 = -rSges(3,3) * t273 + t270 * t298;
t197 = t254 * t316 + t319;
t196 = -t254 * t318 + t317;
t195 = t254 * t317 - t318;
t194 = -t254 * t319 - t316;
t192 = rSges(4,3) * t270 + t273 * t297;
t191 = -rSges(4,3) * t273 + t270 * t297;
t190 = t249 * t270 + t250 * t320;
t189 = -t249 * t320 + t250 * t270;
t188 = -t249 * t273 + t250 * t321;
t187 = -t249 * t321 - t250 * t273;
t177 = (-qJD(6) + t305) * t254 + t255;
t176 = -rSges(5,3) * t254 + (rSges(5,1) * t271 - rSges(5,2) * t268) * t253;
t175 = V_base(5) * rSges(2,3) - t243 * t255 + t304;
t174 = t244 * t255 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t173 = -Icges(5,5) * t254 + (Icges(5,1) * t271 - Icges(5,4) * t268) * t253;
t172 = -Icges(5,6) * t254 + (Icges(5,4) * t271 - Icges(5,2) * t268) * t253;
t171 = -Icges(5,3) * t254 + (Icges(5,5) * t271 - Icges(5,6) * t268) * t253;
t167 = t243 * V_base(4) - t244 * V_base(5) + V_base(3);
t165 = -rSges(6,3) * t254 + (rSges(6,1) * t259 - rSges(6,2) * t258) * t253;
t164 = -Icges(6,5) * t254 + (Icges(6,1) * t259 - Icges(6,4) * t258) * t253;
t163 = -Icges(6,6) * t254 + (Icges(6,4) * t259 - Icges(6,2) * t258) * t253;
t162 = -Icges(6,3) * t254 + (Icges(6,5) * t259 - Icges(6,6) * t258) * t253;
t161 = -rSges(7,3) * t254 + (rSges(7,1) * t250 - rSges(7,2) * t249) * t253;
t160 = -Icges(7,5) * t254 + (Icges(7,1) * t250 - Icges(7,4) * t249) * t253;
t159 = -Icges(7,6) * t254 + (Icges(7,4) * t250 - Icges(7,2) * t249) * t253;
t158 = -Icges(7,3) * t254 + (Icges(7,5) * t250 - Icges(7,6) * t249) * t253;
t155 = t273 * t306 + t170;
t154 = t270 * t306 + t169;
t153 = -pkin(10) * t254 + t253 * t310;
t152 = rSges(5,1) * t214 + rSges(5,2) * t213 + rSges(5,3) * t322;
t151 = rSges(5,1) * t212 + rSges(5,2) * t211 + rSges(5,3) * t323;
t150 = Icges(5,1) * t214 + Icges(5,4) * t213 + Icges(5,5) * t322;
t149 = Icges(5,1) * t212 + Icges(5,4) * t211 + Icges(5,5) * t323;
t148 = Icges(5,4) * t214 + Icges(5,2) * t213 + Icges(5,6) * t322;
t147 = Icges(5,4) * t212 + Icges(5,2) * t211 + Icges(5,6) * t323;
t146 = Icges(5,5) * t214 + Icges(5,6) * t213 + Icges(5,3) * t322;
t145 = Icges(5,5) * t212 + Icges(5,6) * t211 + Icges(5,3) * t323;
t141 = rSges(6,1) * t197 + rSges(6,2) * t196 + rSges(6,3) * t322;
t140 = rSges(6,1) * t195 + rSges(6,2) * t194 + rSges(6,3) * t323;
t139 = Icges(6,1) * t197 + Icges(6,4) * t196 + Icges(6,5) * t322;
t138 = Icges(6,1) * t195 + Icges(6,4) * t194 + Icges(6,5) * t323;
t137 = Icges(6,4) * t197 + Icges(6,2) * t196 + Icges(6,6) * t322;
t136 = Icges(6,4) * t195 + Icges(6,2) * t194 + Icges(6,6) * t323;
t135 = Icges(6,5) * t197 + Icges(6,6) * t196 + Icges(6,3) * t322;
t134 = Icges(6,5) * t195 + Icges(6,6) * t194 + Icges(6,3) * t323;
t133 = rSges(7,1) * t190 + rSges(7,2) * t189 + rSges(7,3) * t322;
t132 = rSges(7,1) * t188 + rSges(7,2) * t187 + rSges(7,3) * t323;
t131 = Icges(7,1) * t190 + Icges(7,4) * t189 + Icges(7,5) * t322;
t130 = Icges(7,1) * t188 + Icges(7,4) * t187 + Icges(7,5) * t323;
t129 = Icges(7,4) * t190 + Icges(7,2) * t189 + Icges(7,6) * t322;
t128 = Icges(7,4) * t188 + Icges(7,2) * t187 + Icges(7,6) * t323;
t127 = Icges(7,5) * t190 + Icges(7,6) * t189 + Icges(7,3) * t322;
t126 = Icges(7,5) * t188 + Icges(7,6) * t187 + Icges(7,3) * t323;
t124 = t242 * t247 + (-t204 - t245) * t255 + t304;
t123 = t205 * t255 - t242 * t248 + t290;
t121 = t270 * t301 + t273 * t284;
t120 = t270 * t284 - t273 * t301;
t119 = t204 * t248 - t205 * t247 + t289;
t118 = t219 * t247 + (-t191 + t311) * t255 + t300;
t117 = t192 * t255 + (-t219 - t334) * t248 + t283;
t116 = t191 * t248 + (-t179 - t192) * t247 + t288;
t115 = -t151 * t224 + t176 * t209 + t282;
t114 = t152 * t224 - t176 * t210 + t278;
t113 = t151 * t210 - t152 * t209 + t281;
t112 = -t140 * t208 + t165 * t169 + t280;
t111 = t141 * t208 - t165 * t170 + t277;
t110 = t140 * t170 - t141 * t169 + t279;
t109 = -t120 * t208 - t132 * t177 + t153 * t169 + t154 * t161 + t280;
t108 = t121 * t208 + t133 * t177 - t153 * t170 - t155 * t161 + t277;
t107 = t120 * t170 - t121 * t169 + t132 * t155 - t133 * t154 + t279;
t1 = t208 * ((-t134 * t169 - t135 * t170 - t162 * t208) * t254 + ((-t137 * t258 + t139 * t259) * t170 + (-t136 * t258 + t138 * t259) * t169 + (-t163 * t258 + t164 * t259) * t208) * t253) / 0.2e1 + t177 * ((-t126 * t154 - t127 * t155 - t158 * t177) * t254 + ((-t129 * t249 + t131 * t250) * t155 + (-t128 * t249 + t130 * t250) * t154 + (-t159 * t249 + t160 * t250) * t177) * t253) / 0.2e1 + m(1) * (t225 ^ 2 + t226 ^ 2 + t227 ^ 2) / 0.2e1 + t210 * ((t146 * t322 + t213 * t148 + t214 * t150) * t210 + (t145 * t322 + t147 * t213 + t149 * t214) * t209 + (t171 * t322 + t172 * t213 + t173 * t214) * t224) / 0.2e1 + t170 * ((t135 * t322 + t196 * t137 + t197 * t139) * t170 + (t134 * t322 + t136 * t196 + t138 * t197) * t169 + (t162 * t322 + t163 * t196 + t164 * t197) * t208) / 0.2e1 + t155 * ((t127 * t322 + t189 * t129 + t190 * t131) * t155 + (t126 * t322 + t128 * t189 + t130 * t190) * t154 + (t158 * t322 + t159 * t189 + t160 * t190) * t177) / 0.2e1 + t209 * ((t146 * t323 + t148 * t211 + t150 * t212) * t210 + (t145 * t323 + t211 * t147 + t212 * t149) * t209 + (t171 * t323 + t172 * t211 + t173 * t212) * t224) / 0.2e1 + t169 * ((t135 * t323 + t137 * t194 + t139 * t195) * t170 + (t134 * t323 + t194 * t136 + t195 * t138) * t169 + (t162 * t323 + t163 * t194 + t164 * t195) * t208) / 0.2e1 + t154 * ((t127 * t323 + t129 * t187 + t131 * t188) * t155 + (t126 * t323 + t187 * t128 + t188 * t130) * t154 + (t158 * t323 + t159 * t187 + t160 * t188) * t177) / 0.2e1 + m(2) * (t167 ^ 2 + t174 ^ 2 + t175 ^ 2) / 0.2e1 + m(3) * (t119 ^ 2 + t123 ^ 2 + t124 ^ 2) / 0.2e1 + m(4) * (t116 ^ 2 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + m(7) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + m(6) * (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + m(5) * (t113 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + t224 * ((-t145 * t209 - t146 * t210 - t171 * t224) * t254 + ((-t148 * t268 + t150 * t271) * t210 + (-t147 * t268 + t149 * t271) * t209 + (-t172 * t268 + t173 * t271) * t224) * t253) / 0.2e1 + (t339 * t270 - t338 * t273) * t247 / 0.2e1 + (t338 * t270 + t339 * t273) * t248 / 0.2e1 + ((-t233 * t270 + t236 * t273 + Icges(1,4)) * V_base(5) + (-t234 * t270 + t237 * t273 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t233 * t273 + t236 * t270 + Icges(1,2)) * V_base(5) + (t234 * t273 + t237 * t270 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t184 * t254 + t186 * t253 + t201 * t272 + t203 * t269) * t248 + (t183 * t254 + t185 * t253 + t200 * t272 + t202 * t269) * t247 + (t217 * t254 + t218 * t253 + t232 * t272 + t235 * t269 + Icges(2,3)) * t255) * t255 / 0.2e1 + t255 * V_base(4) * (Icges(2,5) * t273 - Icges(2,6) * t270) + V_base(5) * t255 * (Icges(2,5) * t270 + Icges(2,6) * t273) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
