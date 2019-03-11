% Calculate kinetic energy for
% S6RRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 15:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPPR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR2_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPPR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR2_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR2_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPPR2_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:24:51
% EndTime: 2019-03-09 15:24:54
% DurationCPUTime: 3.60s
% Computational Cost: add. (1897->312), mult. (1778->438), div. (0->0), fcn. (1560->10), ass. (0->166)
t332 = Icges(5,4) + Icges(6,6);
t331 = Icges(5,1) + Icges(6,2);
t330 = -Icges(5,2) - Icges(6,3);
t235 = qJ(2) + qJ(3);
t224 = pkin(10) + t235;
t222 = cos(t224);
t329 = t332 * t222;
t221 = sin(t224);
t328 = t332 * t221;
t327 = Icges(6,4) - Icges(5,5);
t326 = Icges(6,5) - Icges(5,6);
t325 = t330 * t221 + t329;
t324 = t331 * t222 - t328;
t238 = sin(qJ(1));
t241 = cos(qJ(1));
t323 = t325 * t238 + t326 * t241;
t322 = -t326 * t238 + t325 * t241;
t321 = t324 * t238 + t327 * t241;
t320 = -t327 * t238 + t324 * t241;
t319 = t330 * t222 - t328;
t318 = t331 * t221 + t329;
t317 = Icges(6,1) + Icges(4,3) + Icges(5,3);
t229 = sin(t235);
t230 = cos(t235);
t316 = Icges(4,5) * t230 - Icges(4,6) * t229 + t326 * t221 - t327 * t222;
t301 = Icges(4,4) * t230;
t268 = -Icges(4,2) * t229 + t301;
t157 = -Icges(4,6) * t241 + t238 * t268;
t158 = Icges(4,6) * t238 + t241 * t268;
t302 = Icges(4,4) * t229;
t271 = Icges(4,1) * t230 - t302;
t159 = -Icges(4,5) * t241 + t238 * t271;
t160 = Icges(4,5) * t238 + t241 * t271;
t191 = Icges(4,2) * t230 + t302;
t192 = Icges(4,1) * t229 + t301;
t195 = V_base(5) + (-qJD(2) - qJD(3)) * t241;
t220 = qJD(2) * t238 + V_base(4);
t196 = qJD(3) * t238 + t220;
t225 = V_base(6) + qJD(1);
t315 = (-t191 * t229 + t192 * t230 + t221 * t319 + t222 * t318) * t225 + (-t158 * t229 + t160 * t230 - t221 * t322 + t222 * t320) * t196 + (-t157 * t229 + t159 * t230 - t221 * t323 + t321 * t222) * t195;
t314 = (Icges(4,5) * t229 + Icges(4,6) * t230 - t327 * t221 - t326 * t222) * t225 + (t238 * t317 + t241 * t316) * t196 + (t238 * t316 - t241 * t317) * t195;
t237 = sin(qJ(2));
t310 = pkin(2) * t237;
t309 = pkin(3) * t229;
t308 = pkin(9) * t221;
t240 = cos(qJ(2));
t307 = t240 * pkin(2);
t305 = Icges(2,4) * t238;
t304 = Icges(3,4) * t237;
t303 = Icges(3,4) * t240;
t296 = t222 * t238;
t295 = t222 * t241;
t236 = sin(qJ(6));
t294 = t236 * t238;
t293 = t236 * t241;
t239 = cos(qJ(6));
t292 = t238 * t239;
t291 = t239 * t241;
t288 = pkin(3) * t230;
t122 = qJ(4) * t238 + t241 * t288;
t273 = pkin(4) * t222 + qJ(5) * t221;
t164 = t273 * t241;
t290 = -t122 - t164;
t151 = -pkin(8) * t241 + t238 * t307;
t217 = t238 * pkin(1) - t241 * pkin(7);
t289 = -t151 - t217;
t286 = qJD(5) * t221;
t285 = qJD(6) * t222;
t284 = V_base(5) * pkin(6) + V_base(1);
t121 = -qJ(4) * t241 + t238 * t288;
t281 = -t121 + t289;
t184 = pkin(4) * t221 - qJ(5) * t222;
t280 = -t184 - t309;
t219 = -qJD(2) * t241 + V_base(5);
t279 = t219 * t310 + t284;
t163 = t273 * t238;
t278 = -t163 + t281;
t277 = rSges(3,1) * t240 - rSges(3,2) * t237;
t276 = rSges(4,1) * t230 - rSges(4,2) * t229;
t275 = rSges(5,1) * t222 - rSges(5,2) * t221;
t274 = -rSges(6,2) * t222 + rSges(6,3) * t221;
t272 = Icges(3,1) * t240 - t304;
t269 = -Icges(3,2) * t237 + t303;
t265 = Icges(3,5) * t240 - Icges(3,6) * t237;
t260 = qJD(4) * t238 + t195 * t309 + t279;
t218 = t241 * pkin(1) + t238 * pkin(7);
t259 = -V_base(4) * pkin(6) + t225 * t218 + V_base(2);
t258 = V_base(4) * t217 - t218 * V_base(5) + V_base(3);
t257 = t195 * t184 + t241 * t286 + t260;
t253 = (-Icges(3,3) * t241 + t238 * t265) * t219 + (Icges(3,3) * t238 + t241 * t265) * t220 + (Icges(3,5) * t237 + Icges(3,6) * t240) * t225;
t152 = pkin(8) * t238 + t241 * t307;
t252 = t220 * t151 - t152 * t219 + t258;
t251 = t225 * t152 - t220 * t310 + t259;
t250 = t196 * t121 + t252;
t249 = -qJD(4) * t241 + t225 * t122 + t251;
t248 = -qJD(5) * t222 + t196 * t163 + t250;
t247 = t225 * t164 + t238 * t286 + t249;
t167 = -Icges(3,6) * t241 + t238 * t269;
t168 = Icges(3,6) * t238 + t241 * t269;
t169 = -Icges(3,5) * t241 + t238 * t272;
t170 = Icges(3,5) * t238 + t241 * t272;
t208 = Icges(3,2) * t240 + t304;
t211 = Icges(3,1) * t237 + t303;
t243 = (-t168 * t237 + t170 * t240) * t220 + (-t167 * t237 + t169 * t240) * t219 + (-t208 * t237 + t211 * t240) * t225;
t231 = Icges(2,4) * t241;
t216 = rSges(2,1) * t241 - rSges(2,2) * t238;
t215 = rSges(2,1) * t238 + rSges(2,2) * t241;
t214 = rSges(3,1) * t237 + rSges(3,2) * t240;
t213 = Icges(2,1) * t241 - t305;
t212 = Icges(2,1) * t238 + t231;
t210 = -Icges(2,2) * t238 + t231;
t209 = Icges(2,2) * t241 + t305;
t202 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t201 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t200 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t194 = qJD(6) * t221 + t225;
t193 = rSges(4,1) * t229 + rSges(4,2) * t230;
t188 = -pkin(5) * t241 + pkin(9) * t296;
t187 = pkin(5) * t238 + pkin(9) * t295;
t186 = rSges(5,1) * t221 + rSges(5,2) * t222;
t185 = -rSges(6,2) * t221 - rSges(6,3) * t222;
t177 = t221 * t294 - t291;
t176 = t221 * t292 + t293;
t175 = t221 * t293 + t292;
t174 = t221 * t291 - t294;
t172 = rSges(3,3) * t238 + t241 * t277;
t171 = -rSges(3,3) * t241 + t238 * t277;
t162 = rSges(4,3) * t238 + t241 * t276;
t161 = -rSges(4,3) * t241 + t238 * t276;
t154 = t241 * t285 + t196;
t153 = t238 * t285 + t195;
t150 = -rSges(6,1) * t241 + t238 * t274;
t149 = rSges(6,1) * t238 + t241 * t274;
t148 = rSges(5,3) * t238 + t241 * t275;
t147 = -rSges(5,3) * t241 + t238 * t275;
t146 = V_base(5) * rSges(2,3) - t215 * t225 + t284;
t145 = t216 * t225 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t131 = t215 * V_base(4) - t216 * V_base(5) + V_base(3);
t129 = rSges(7,3) * t221 + (-rSges(7,1) * t236 - rSges(7,2) * t239) * t222;
t127 = Icges(7,5) * t221 + (-Icges(7,1) * t236 - Icges(7,4) * t239) * t222;
t126 = Icges(7,6) * t221 + (-Icges(7,4) * t236 - Icges(7,2) * t239) * t222;
t125 = Icges(7,3) * t221 + (-Icges(7,5) * t236 - Icges(7,6) * t239) * t222;
t118 = rSges(7,1) * t177 + rSges(7,2) * t176 + rSges(7,3) * t296;
t117 = rSges(7,1) * t175 + rSges(7,2) * t174 + rSges(7,3) * t295;
t116 = Icges(7,1) * t177 + Icges(7,4) * t176 + Icges(7,5) * t296;
t115 = Icges(7,1) * t175 + Icges(7,4) * t174 + Icges(7,5) * t295;
t114 = Icges(7,4) * t177 + Icges(7,2) * t176 + Icges(7,6) * t296;
t113 = Icges(7,4) * t175 + Icges(7,2) * t174 + Icges(7,6) * t295;
t112 = Icges(7,5) * t177 + Icges(7,6) * t176 + Icges(7,3) * t296;
t111 = Icges(7,5) * t175 + Icges(7,6) * t174 + Icges(7,3) * t295;
t110 = t214 * t219 + (-t171 - t217) * t225 + t284;
t109 = t172 * t225 - t214 * t220 + t259;
t108 = t171 * t220 - t172 * t219 + t258;
t107 = t193 * t195 + (-t161 + t289) * t225 + t279;
t106 = t162 * t225 - t193 * t196 + t251;
t105 = t161 * t196 - t162 * t195 + t252;
t104 = t186 * t195 + (-t147 + t281) * t225 + t260;
t103 = t148 * t225 + (-t186 - t309) * t196 + t249;
t102 = t185 * t195 + (-t150 + t278) * t225 + t257;
t101 = t149 * t225 + (-t185 + t280) * t196 + t247;
t100 = t147 * t196 + (-t122 - t148) * t195 + t250;
t99 = t150 * t196 + (-t149 + t290) * t195 + t248;
t98 = t195 * t308 - t118 * t194 + t129 * t153 + (-t188 + t278) * t225 + t257;
t97 = t117 * t194 - t129 * t154 + t187 * t225 + (t280 - t308) * t196 + t247;
t96 = -t117 * t153 + t118 * t154 + t188 * t196 + (-t187 + t290) * t195 + t248;
t1 = m(1) * (t200 ^ 2 + t201 ^ 2 + t202 ^ 2) / 0.2e1 + m(2) * (t131 ^ 2 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + m(4) * (t105 ^ 2 + t106 ^ 2 + t107 ^ 2) / 0.2e1 + m(3) * (t108 ^ 2 + t109 ^ 2 + t110 ^ 2) / 0.2e1 + m(7) * (t96 ^ 2 + t97 ^ 2 + t98 ^ 2) / 0.2e1 + m(6) * (t101 ^ 2 + t102 ^ 2 + t99 ^ 2) / 0.2e1 + m(5) * (t100 ^ 2 + t103 ^ 2 + t104 ^ 2) / 0.2e1 + t194 * ((t111 * t154 + t112 * t153 + t125 * t194) * t221 + ((-t113 * t239 - t115 * t236) * t154 + (-t114 * t239 - t116 * t236) * t153 + (-t126 * t239 - t127 * t236) * t194) * t222) / 0.2e1 + t154 * ((t111 * t295 + t174 * t113 + t175 * t115) * t154 + (t112 * t295 + t114 * t174 + t116 * t175) * t153 + (t125 * t295 + t126 * t174 + t127 * t175) * t194) / 0.2e1 + t153 * ((t111 * t296 + t113 * t176 + t115 * t177) * t154 + (t112 * t296 + t176 * t114 + t177 * t116) * t153 + (t125 * t296 + t126 * t176 + t127 * t177) * t194) / 0.2e1 + t220 * (t253 * t238 + t243 * t241) / 0.2e1 + t219 * (t243 * t238 - t253 * t241) / 0.2e1 + ((-t209 * t238 + t212 * t241 + Icges(1,4)) * V_base(5) + (-t210 * t238 + t213 * t241 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t209 * t241 + t212 * t238 + Icges(1,2)) * V_base(5) + (t210 * t241 + t213 * t238 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (t315 * t238 - t314 * t241) * t195 / 0.2e1 + (t314 * t238 + t315 * t241) * t196 / 0.2e1 + ((t168 * t240 + t170 * t237) * t220 + (t167 * t240 + t169 * t237) * t219 + (t158 * t230 + t160 * t229 + t221 * t320 + t222 * t322) * t196 + (t157 * t230 + t159 * t229 + t321 * t221 + t222 * t323) * t195 + (t191 * t230 + t192 * t229 + t208 * t240 + t211 * t237 + t221 * t318 - t319 * t222 + Icges(2,3)) * t225) * t225 / 0.2e1 + t225 * V_base(4) * (Icges(2,5) * t241 - Icges(2,6) * t238) + t225 * V_base(5) * (Icges(2,5) * t238 + Icges(2,6) * t241) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
