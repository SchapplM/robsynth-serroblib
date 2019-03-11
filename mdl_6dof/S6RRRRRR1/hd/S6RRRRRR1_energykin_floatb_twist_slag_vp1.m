% Calculate kinetic energy for
% S6RRRRRR1
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
% Datum: 2019-03-10 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR1_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR1_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRRR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR1_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR1_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR1_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRR1_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:27:47
% EndTime: 2019-03-10 03:27:49
% DurationCPUTime: 2.36s
% Computational Cost: add. (2281->338), mult. (1880->511), div. (0->0), fcn. (1662->12), ass. (0->183)
t248 = sin(qJ(2));
t322 = pkin(2) * t248;
t246 = qJ(2) + qJ(3);
t238 = sin(t246);
t321 = pkin(3) * t238;
t241 = qJ(4) + t246;
t229 = sin(t241);
t320 = pkin(4) * t229;
t251 = cos(qJ(2));
t319 = t251 * pkin(2);
t249 = sin(qJ(1));
t317 = Icges(2,4) * t249;
t316 = Icges(3,4) * t248;
t315 = Icges(3,4) * t251;
t314 = Icges(4,4) * t238;
t239 = cos(t246);
t313 = Icges(4,4) * t239;
t312 = Icges(5,4) * t229;
t230 = cos(t241);
t311 = Icges(5,4) * t230;
t232 = qJ(5) + t241;
t226 = sin(t232);
t310 = Icges(6,4) * t226;
t227 = cos(t232);
t309 = Icges(6,4) * t227;
t308 = t226 * t249;
t252 = cos(qJ(1));
t307 = t226 * t252;
t247 = sin(qJ(6));
t306 = t247 * t249;
t305 = t247 * t252;
t250 = cos(qJ(6));
t304 = t249 * t250;
t303 = t250 * t252;
t157 = -pkin(8) * t252 + t249 * t319;
t222 = t249 * pkin(1) - t252 * pkin(7);
t302 = -t157 - t222;
t301 = pkin(4) * t230;
t300 = pkin(3) * t239;
t297 = qJD(6) * t226;
t296 = -qJD(2) - qJD(3);
t295 = V_base(5) * pkin(6) + V_base(1);
t127 = -pkin(9) * t252 + t249 * t300;
t292 = -t127 + t302;
t225 = qJD(2) * t249 + V_base(4);
t233 = V_base(6) + qJD(1);
t291 = -qJD(4) + t296;
t224 = -qJD(2) * t252 + V_base(5);
t290 = t224 * t322 + t295;
t123 = -pkin(10) * t252 + t249 * t301;
t289 = -t123 + t292;
t203 = qJD(3) * t249 + t225;
t202 = t252 * t296 + V_base(5);
t288 = t202 * t321 + t290;
t287 = pkin(5) * t227 + pkin(11) * t226;
t286 = rSges(3,1) * t251 - rSges(3,2) * t248;
t285 = rSges(4,1) * t239 - rSges(4,2) * t238;
t284 = rSges(5,1) * t230 - rSges(5,2) * t229;
t283 = rSges(6,1) * t227 - rSges(6,2) * t226;
t191 = qJD(4) * t249 + t203;
t282 = Icges(3,1) * t251 - t316;
t281 = Icges(4,1) * t239 - t314;
t280 = Icges(5,1) * t230 - t312;
t279 = Icges(6,1) * t227 - t310;
t278 = -Icges(3,2) * t248 + t315;
t277 = -Icges(4,2) * t238 + t313;
t276 = -Icges(5,2) * t229 + t311;
t275 = -Icges(6,2) * t226 + t309;
t274 = Icges(3,5) * t251 - Icges(3,6) * t248;
t273 = Icges(4,5) * t239 - Icges(4,6) * t238;
t272 = Icges(5,5) * t230 - Icges(5,6) * t229;
t271 = Icges(6,5) * t227 - Icges(6,6) * t226;
t190 = t252 * t291 + V_base(5);
t270 = t190 * t320 + t288;
t223 = t252 * pkin(1) + t249 * pkin(7);
t269 = -V_base(4) * pkin(6) + t233 * t223 + V_base(2);
t170 = qJD(5) * t249 + t191;
t268 = V_base(4) * t222 - t223 * V_base(5) + V_base(3);
t169 = V_base(5) + (-qJD(5) + t291) * t252;
t267 = (-Icges(6,3) * t252 + t249 * t271) * t169 + (Icges(6,3) * t249 + t252 * t271) * t170 + (Icges(6,5) * t226 + Icges(6,6) * t227) * t233;
t266 = (-Icges(5,3) * t252 + t249 * t272) * t190 + (Icges(5,3) * t249 + t252 * t272) * t191 + (Icges(5,5) * t229 + Icges(5,6) * t230) * t233;
t265 = (-Icges(4,3) * t252 + t249 * t273) * t202 + (Icges(4,3) * t249 + t252 * t273) * t203 + (Icges(4,5) * t238 + Icges(4,6) * t239) * t233;
t264 = (-Icges(3,3) * t252 + t249 * t274) * t224 + (Icges(3,3) * t249 + t252 * t274) * t225 + (Icges(3,5) * t248 + Icges(3,6) * t251) * t233;
t158 = pkin(8) * t249 + t252 * t319;
t263 = t225 * t157 - t158 * t224 + t268;
t262 = t233 * t158 - t225 * t322 + t269;
t128 = pkin(9) * t249 + t252 * t300;
t261 = t203 * t127 - t128 * t202 + t263;
t260 = t233 * t128 - t203 * t321 + t262;
t124 = pkin(10) * t249 + t252 * t301;
t259 = t191 * t123 - t124 * t190 + t261;
t258 = t233 * t124 - t191 * t320 + t260;
t139 = -Icges(6,6) * t252 + t249 * t275;
t140 = Icges(6,6) * t249 + t252 * t275;
t141 = -Icges(6,5) * t252 + t249 * t279;
t142 = Icges(6,5) * t249 + t252 * t279;
t186 = Icges(6,2) * t227 + t310;
t187 = Icges(6,1) * t226 + t309;
t257 = (-t140 * t226 + t142 * t227) * t170 + (-t139 * t226 + t141 * t227) * t169 + (-t186 * t226 + t187 * t227) * t233;
t150 = -Icges(5,6) * t252 + t249 * t276;
t151 = Icges(5,6) * t249 + t252 * t276;
t152 = -Icges(5,5) * t252 + t249 * t280;
t153 = Icges(5,5) * t249 + t252 * t280;
t193 = Icges(5,2) * t230 + t312;
t194 = Icges(5,1) * t229 + t311;
t256 = (-t151 * t229 + t153 * t230) * t191 + (-t150 * t229 + t152 * t230) * t190 + (-t193 * t229 + t194 * t230) * t233;
t161 = -Icges(4,6) * t252 + t249 * t277;
t162 = Icges(4,6) * t249 + t252 * t277;
t163 = -Icges(4,5) * t252 + t249 * t281;
t164 = Icges(4,5) * t249 + t252 * t281;
t199 = Icges(4,2) * t239 + t314;
t200 = Icges(4,1) * t238 + t313;
t255 = (-t162 * t238 + t164 * t239) * t203 + (-t161 * t238 + t163 * t239) * t202 + (-t199 * t238 + t200 * t239) * t233;
t173 = -Icges(3,6) * t252 + t249 * t278;
t174 = Icges(3,6) * t249 + t252 * t278;
t175 = -Icges(3,5) * t252 + t249 * t282;
t176 = Icges(3,5) * t249 + t252 * t282;
t213 = Icges(3,2) * t251 + t316;
t216 = Icges(3,1) * t248 + t315;
t254 = (-t174 * t248 + t176 * t251) * t225 + (-t173 * t248 + t175 * t251) * t224 + (-t213 * t248 + t216 * t251) * t233;
t240 = Icges(2,4) * t252;
t221 = rSges(2,1) * t252 - rSges(2,2) * t249;
t220 = rSges(2,1) * t249 + rSges(2,2) * t252;
t219 = rSges(3,1) * t248 + rSges(3,2) * t251;
t218 = Icges(2,1) * t252 - t317;
t217 = Icges(2,1) * t249 + t240;
t215 = -Icges(2,2) * t249 + t240;
t214 = Icges(2,2) * t252 + t317;
t209 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t208 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t207 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t201 = rSges(4,1) * t238 + rSges(4,2) * t239;
t197 = -qJD(6) * t227 + t233;
t195 = rSges(5,1) * t229 + rSges(5,2) * t230;
t189 = pkin(5) * t226 - pkin(11) * t227;
t188 = rSges(6,1) * t226 + rSges(6,2) * t227;
t183 = rSges(3,3) * t249 + t252 * t286;
t182 = -rSges(3,3) * t252 + t249 * t286;
t180 = t227 * t303 + t306;
t179 = -t227 * t305 + t304;
t178 = t227 * t304 - t305;
t177 = -t227 * t306 - t303;
t168 = rSges(4,3) * t249 + t252 * t285;
t167 = -rSges(4,3) * t252 + t249 * t285;
t166 = t287 * t252;
t165 = t287 * t249;
t155 = rSges(5,3) * t249 + t252 * t284;
t154 = -rSges(5,3) * t252 + t249 * t284;
t147 = V_base(5) * rSges(2,3) - t220 * t233 + t295;
t146 = t221 * t233 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t145 = rSges(6,3) * t249 + t252 * t283;
t144 = -rSges(6,3) * t252 + t249 * t283;
t143 = t220 * V_base(4) - t221 * V_base(5) + V_base(3);
t135 = -rSges(7,3) * t227 + (rSges(7,1) * t250 - rSges(7,2) * t247) * t226;
t133 = -Icges(7,5) * t227 + (Icges(7,1) * t250 - Icges(7,4) * t247) * t226;
t132 = -Icges(7,6) * t227 + (Icges(7,4) * t250 - Icges(7,2) * t247) * t226;
t131 = -Icges(7,3) * t227 + (Icges(7,5) * t250 - Icges(7,6) * t247) * t226;
t130 = t252 * t297 + t170;
t129 = t249 * t297 + t169;
t122 = rSges(7,1) * t180 + rSges(7,2) * t179 + rSges(7,3) * t307;
t121 = rSges(7,1) * t178 + rSges(7,2) * t177 + rSges(7,3) * t308;
t120 = Icges(7,1) * t180 + Icges(7,4) * t179 + Icges(7,5) * t307;
t119 = Icges(7,1) * t178 + Icges(7,4) * t177 + Icges(7,5) * t308;
t118 = Icges(7,4) * t180 + Icges(7,2) * t179 + Icges(7,6) * t307;
t117 = Icges(7,4) * t178 + Icges(7,2) * t177 + Icges(7,6) * t308;
t116 = Icges(7,5) * t180 + Icges(7,6) * t179 + Icges(7,3) * t307;
t115 = Icges(7,5) * t178 + Icges(7,6) * t177 + Icges(7,3) * t308;
t113 = t219 * t224 + (-t182 - t222) * t233 + t295;
t112 = t183 * t233 - t219 * t225 + t269;
t110 = t182 * t225 - t183 * t224 + t268;
t109 = t201 * t202 + (-t167 + t302) * t233 + t290;
t108 = t168 * t233 - t201 * t203 + t262;
t107 = t167 * t203 - t168 * t202 + t263;
t106 = t190 * t195 + (-t154 + t292) * t233 + t288;
t105 = t155 * t233 - t191 * t195 + t260;
t104 = t154 * t191 - t155 * t190 + t261;
t103 = t169 * t188 + (-t144 + t289) * t233 + t270;
t102 = t145 * t233 - t170 * t188 + t258;
t101 = -t121 * t197 + t129 * t135 + t169 * t189 + (-t165 + t289) * t233 + t270;
t100 = t122 * t197 - t130 * t135 + t166 * t233 - t170 * t189 + t258;
t99 = t144 * t170 - t145 * t169 + t259;
t98 = t121 * t130 - t122 * t129 + t165 * t170 - t166 * t169 + t259;
t1 = m(1) * (t207 ^ 2 + t208 ^ 2 + t209 ^ 2) / 0.2e1 + m(2) * (t143 ^ 2 + t146 ^ 2 + t147 ^ 2) / 0.2e1 + t130 * ((t116 * t307 + t179 * t118 + t180 * t120) * t130 + (t115 * t307 + t117 * t179 + t119 * t180) * t129 + (t131 * t307 + t132 * t179 + t133 * t180) * t197) / 0.2e1 + t129 * ((t116 * t308 + t118 * t177 + t120 * t178) * t130 + (t115 * t308 + t177 * t117 + t178 * t119) * t129 + (t131 * t308 + t132 * t177 + t133 * t178) * t197) / 0.2e1 + t225 * (t264 * t249 + t254 * t252) / 0.2e1 + t224 * (t254 * t249 - t264 * t252) / 0.2e1 + t203 * (t265 * t249 + t255 * t252) / 0.2e1 + t202 * (t255 * t249 - t265 * t252) / 0.2e1 + t191 * (t266 * t249 + t256 * t252) / 0.2e1 + t190 * (t256 * t249 - t266 * t252) / 0.2e1 + t170 * (t267 * t249 + t257 * t252) / 0.2e1 + t169 * (t257 * t249 - t267 * t252) / 0.2e1 + t197 * ((-t115 * t129 - t116 * t130 - t131 * t197) * t227 + ((-t118 * t247 + t120 * t250) * t130 + (-t117 * t247 + t119 * t250) * t129 + (-t132 * t247 + t133 * t250) * t197) * t226) / 0.2e1 + m(3) * (t110 ^ 2 + t112 ^ 2 + t113 ^ 2) / 0.2e1 + m(7) * (t100 ^ 2 + t101 ^ 2 + t98 ^ 2) / 0.2e1 + m(6) * (t102 ^ 2 + t103 ^ 2 + t99 ^ 2) / 0.2e1 + m(5) * (t104 ^ 2 + t105 ^ 2 + t106 ^ 2) / 0.2e1 + m(4) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + ((-t214 * t249 + t217 * t252 + Icges(1,4)) * V_base(5) + (-t249 * t215 + t252 * t218 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t252 * t214 + t249 * t217 + Icges(1,2)) * V_base(5) + (t215 * t252 + t218 * t249 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t140 * t227 + t142 * t226) * t170 + (t139 * t227 + t141 * t226) * t169 + (t174 * t251 + t176 * t248) * t225 + (t173 * t251 + t175 * t248) * t224 + (t162 * t239 + t164 * t238) * t203 + (t161 * t239 + t163 * t238) * t202 + (t151 * t230 + t153 * t229) * t191 + (t150 * t230 + t152 * t229) * t190 + (t227 * t186 + t226 * t187 + t230 * t193 + t229 * t194 + t239 * t199 + t238 * t200 + t251 * t213 + t248 * t216 + Icges(2,3)) * t233) * t233 / 0.2e1 + V_base(4) * t233 * (Icges(2,5) * t252 - Icges(2,6) * t249) + V_base(5) * t233 * (Icges(2,5) * t249 + Icges(2,6) * t252) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
