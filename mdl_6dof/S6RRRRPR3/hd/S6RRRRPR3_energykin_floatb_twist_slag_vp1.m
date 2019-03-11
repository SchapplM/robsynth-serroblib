% Calculate kinetic energy for
% S6RRRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-09 22:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR3_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR3_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRPR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR3_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR3_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR3_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPR3_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:01:49
% EndTime: 2019-03-09 22:01:53
% DurationCPUTime: 3.74s
% Computational Cost: add. (1939->316), mult. (1820->460), div. (0->0), fcn. (1602->10), ass. (0->170)
t333 = Icges(5,4) + Icges(6,6);
t332 = Icges(5,1) + Icges(6,2);
t331 = -Icges(5,2) - Icges(6,3);
t237 = qJ(2) + qJ(3);
t233 = qJ(4) + t237;
t224 = cos(t233);
t330 = t333 * t224;
t223 = sin(t233);
t329 = t333 * t223;
t328 = Icges(6,4) - Icges(5,5);
t327 = Icges(6,5) - Icges(5,6);
t326 = t331 * t223 + t330;
t325 = t332 * t224 - t329;
t324 = Icges(6,1) + Icges(5,3);
t240 = sin(qJ(1));
t243 = cos(qJ(1));
t323 = t326 * t240 + t327 * t243;
t322 = -t327 * t240 + t326 * t243;
t321 = t325 * t240 + t328 * t243;
t320 = -t328 * t240 + t325 * t243;
t319 = t331 * t224 - t329;
t318 = t332 * t223 + t330;
t317 = t327 * t223 - t328 * t224;
t286 = -qJD(2) - qJD(3);
t178 = V_base(5) + (-qJD(4) + t286) * t243;
t222 = qJD(2) * t240 + V_base(4);
t197 = qJD(3) * t240 + t222;
t179 = qJD(4) * t240 + t197;
t226 = V_base(6) + qJD(1);
t316 = (t319 * t223 + t318 * t224) * t226 + (-t322 * t223 + t320 * t224) * t179 + (-t323 * t223 + t321 * t224) * t178;
t315 = (-t328 * t223 - t327 * t224) * t226 + (t324 * t240 + t317 * t243) * t179 + (t317 * t240 - t324 * t243) * t178;
t239 = sin(qJ(2));
t311 = pkin(2) * t239;
t230 = sin(t237);
t310 = pkin(3) * t230;
t309 = pkin(10) * t223;
t242 = cos(qJ(2));
t308 = t242 * pkin(2);
t306 = Icges(2,4) * t240;
t305 = Icges(3,4) * t239;
t304 = Icges(3,4) * t242;
t303 = Icges(4,4) * t230;
t231 = cos(t237);
t302 = Icges(4,4) * t231;
t297 = t224 * t240;
t296 = t224 * t243;
t238 = sin(qJ(6));
t295 = t238 * t240;
t294 = t238 * t243;
t241 = cos(qJ(6));
t293 = t240 * t241;
t292 = t241 * t243;
t153 = -pkin(8) * t243 + t240 * t308;
t219 = t240 * pkin(1) - t243 * pkin(7);
t291 = -t153 - t219;
t290 = pkin(3) * t231;
t288 = qJD(5) * t223;
t287 = qJD(6) * t224;
t285 = V_base(5) * pkin(6) + V_base(1);
t122 = -pkin(9) * t243 + t240 * t290;
t282 = -t122 + t291;
t221 = -qJD(2) * t243 + V_base(5);
t281 = t221 * t311 + t285;
t274 = pkin(4) * t224 + qJ(5) * t223;
t163 = t274 * t240;
t280 = -t163 + t282;
t196 = t243 * t286 + V_base(5);
t279 = t196 * t310 + t281;
t278 = rSges(3,1) * t242 - rSges(3,2) * t239;
t277 = rSges(4,1) * t231 - rSges(4,2) * t230;
t276 = rSges(5,1) * t224 - rSges(5,2) * t223;
t275 = -rSges(6,2) * t224 + rSges(6,3) * t223;
t273 = Icges(3,1) * t242 - t305;
t272 = Icges(4,1) * t231 - t303;
t270 = -Icges(3,2) * t239 + t304;
t269 = -Icges(4,2) * t230 + t302;
t266 = Icges(3,5) * t242 - Icges(3,6) * t239;
t265 = Icges(4,5) * t231 - Icges(4,6) * t230;
t220 = t243 * pkin(1) + t240 * pkin(7);
t261 = -V_base(4) * pkin(6) + t226 * t220 + V_base(2);
t260 = V_base(4) * t219 - t220 * V_base(5) + V_base(3);
t186 = pkin(4) * t223 - qJ(5) * t224;
t259 = t178 * t186 + t243 * t288 + t279;
t256 = (-Icges(4,3) * t243 + t240 * t265) * t196 + (Icges(4,3) * t240 + t243 * t265) * t197 + (Icges(4,5) * t230 + Icges(4,6) * t231) * t226;
t255 = (-Icges(3,3) * t243 + t240 * t266) * t221 + (Icges(3,3) * t240 + t243 * t266) * t222 + (Icges(3,5) * t239 + Icges(3,6) * t242) * t226;
t154 = pkin(8) * t240 + t243 * t308;
t254 = t222 * t153 - t154 * t221 + t260;
t253 = t226 * t154 - t222 * t311 + t261;
t123 = pkin(9) * t240 + t243 * t290;
t252 = t197 * t122 - t123 * t196 + t254;
t251 = t226 * t123 - t197 * t310 + t253;
t164 = t274 * t243;
t250 = t226 * t164 + t240 * t288 + t251;
t249 = -qJD(5) * t224 + t179 * t163 + t252;
t157 = -Icges(4,6) * t243 + t240 * t269;
t158 = Icges(4,6) * t240 + t243 * t269;
t159 = -Icges(4,5) * t243 + t240 * t272;
t160 = Icges(4,5) * t240 + t243 * t272;
t193 = Icges(4,2) * t231 + t303;
t194 = Icges(4,1) * t230 + t302;
t246 = (-t158 * t230 + t160 * t231) * t197 + (-t157 * t230 + t159 * t231) * t196 + (-t193 * t230 + t194 * t231) * t226;
t167 = -Icges(3,6) * t243 + t240 * t270;
t168 = Icges(3,6) * t240 + t243 * t270;
t169 = -Icges(3,5) * t243 + t240 * t273;
t170 = Icges(3,5) * t240 + t243 * t273;
t210 = Icges(3,2) * t242 + t305;
t213 = Icges(3,1) * t239 + t304;
t245 = (-t168 * t239 + t170 * t242) * t222 + (-t167 * t239 + t169 * t242) * t221 + (-t210 * t239 + t213 * t242) * t226;
t232 = Icges(2,4) * t243;
t218 = rSges(2,1) * t243 - rSges(2,2) * t240;
t217 = rSges(2,1) * t240 + rSges(2,2) * t243;
t216 = rSges(3,1) * t239 + rSges(3,2) * t242;
t215 = Icges(2,1) * t243 - t306;
t214 = Icges(2,1) * t240 + t232;
t212 = -Icges(2,2) * t240 + t232;
t211 = Icges(2,2) * t243 + t306;
t204 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t203 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t202 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t198 = qJD(6) * t223 + t226;
t195 = rSges(4,1) * t230 + rSges(4,2) * t231;
t190 = -pkin(5) * t243 + pkin(10) * t297;
t189 = pkin(5) * t240 + pkin(10) * t296;
t188 = rSges(5,1) * t223 + rSges(5,2) * t224;
t187 = -rSges(6,2) * t223 - rSges(6,3) * t224;
t177 = t223 * t295 - t292;
t176 = t223 * t293 + t294;
t175 = t223 * t294 + t293;
t174 = t223 * t292 - t295;
t172 = rSges(3,3) * t240 + t243 * t278;
t171 = -rSges(3,3) * t243 + t240 * t278;
t162 = rSges(4,3) * t240 + t243 * t277;
t161 = -rSges(4,3) * t243 + t240 * t277;
t152 = -rSges(6,1) * t243 + t240 * t275;
t151 = rSges(6,1) * t240 + t243 * t275;
t150 = rSges(5,3) * t240 + t243 * t276;
t149 = -rSges(5,3) * t243 + t240 * t276;
t136 = V_base(5) * rSges(2,3) - t217 * t226 + t285;
t135 = t218 * t226 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t133 = t243 * t287 + t179;
t132 = t240 * t287 + t178;
t131 = t217 * V_base(4) - t218 * V_base(5) + V_base(3);
t130 = rSges(7,3) * t223 + (-rSges(7,1) * t238 - rSges(7,2) * t241) * t224;
t129 = Icges(7,5) * t223 + (-Icges(7,1) * t238 - Icges(7,4) * t241) * t224;
t128 = Icges(7,6) * t223 + (-Icges(7,4) * t238 - Icges(7,2) * t241) * t224;
t127 = Icges(7,3) * t223 + (-Icges(7,5) * t238 - Icges(7,6) * t241) * t224;
t118 = rSges(7,1) * t177 + rSges(7,2) * t176 + rSges(7,3) * t297;
t117 = rSges(7,1) * t175 + rSges(7,2) * t174 + rSges(7,3) * t296;
t116 = Icges(7,1) * t177 + Icges(7,4) * t176 + Icges(7,5) * t297;
t115 = Icges(7,1) * t175 + Icges(7,4) * t174 + Icges(7,5) * t296;
t114 = Icges(7,4) * t177 + Icges(7,2) * t176 + Icges(7,6) * t297;
t113 = Icges(7,4) * t175 + Icges(7,2) * t174 + Icges(7,6) * t296;
t112 = Icges(7,5) * t177 + Icges(7,6) * t176 + Icges(7,3) * t297;
t111 = Icges(7,5) * t175 + Icges(7,6) * t174 + Icges(7,3) * t296;
t110 = t216 * t221 + (-t171 - t219) * t226 + t285;
t109 = t172 * t226 - t216 * t222 + t261;
t108 = t171 * t222 - t172 * t221 + t260;
t107 = t195 * t196 + (-t161 + t291) * t226 + t281;
t106 = t162 * t226 - t195 * t197 + t253;
t105 = t161 * t197 - t162 * t196 + t254;
t104 = t178 * t188 + (-t149 + t282) * t226 + t279;
t103 = t150 * t226 - t179 * t188 + t251;
t102 = t178 * t187 + (-t152 + t280) * t226 + t259;
t101 = t151 * t226 + (-t186 - t187) * t179 + t250;
t100 = t149 * t179 - t150 * t178 + t252;
t99 = t178 * t309 - t118 * t198 + t130 * t132 + (-t190 + t280) * t226 + t259;
t98 = t117 * t198 - t130 * t133 + t189 * t226 + (-t186 - t309) * t179 + t250;
t97 = t152 * t179 + (-t151 - t164) * t178 + t249;
t96 = -t117 * t132 + t118 * t133 + t179 * t190 + (-t164 - t189) * t178 + t249;
t1 = t198 * ((t111 * t133 + t112 * t132 + t127 * t198) * t223 + ((-t113 * t241 - t115 * t238) * t133 + (-t114 * t241 - t116 * t238) * t132 + (-t128 * t241 - t129 * t238) * t198) * t224) / 0.2e1 + m(2) * (t131 ^ 2 + t135 ^ 2 + t136 ^ 2) / 0.2e1 + m(3) * (t108 ^ 2 + t109 ^ 2 + t110 ^ 2) / 0.2e1 + m(6) * (t101 ^ 2 + t102 ^ 2 + t97 ^ 2) / 0.2e1 + m(5) * (t100 ^ 2 + t103 ^ 2 + t104 ^ 2) / 0.2e1 + m(4) * (t105 ^ 2 + t106 ^ 2 + t107 ^ 2) / 0.2e1 + m(7) * (t96 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(1) * (t202 ^ 2 + t203 ^ 2 + t204 ^ 2) / 0.2e1 + t222 * (t255 * t240 + t245 * t243) / 0.2e1 + t221 * (t245 * t240 - t255 * t243) / 0.2e1 + t197 * (t256 * t240 + t246 * t243) / 0.2e1 + t196 * (t246 * t240 - t256 * t243) / 0.2e1 + t133 * ((t111 * t296 + t174 * t113 + t175 * t115) * t133 + (t112 * t296 + t114 * t174 + t116 * t175) * t132 + (t127 * t296 + t128 * t174 + t129 * t175) * t198) / 0.2e1 + t132 * ((t111 * t297 + t113 * t176 + t115 * t177) * t133 + (t112 * t297 + t176 * t114 + t177 * t116) * t132 + (t127 * t297 + t128 * t176 + t129 * t177) * t198) / 0.2e1 + (t316 * t240 - t315 * t243) * t178 / 0.2e1 + (t315 * t240 + t316 * t243) * t179 / 0.2e1 + ((-t211 * t240 + t214 * t243 + Icges(1,4)) * V_base(5) + (-t212 * t240 + t215 * t243 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t211 * t243 + t214 * t240 + Icges(1,2)) * V_base(5) + (t212 * t243 + t215 * t240 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t168 * t242 + t170 * t239) * t222 + (t167 * t242 + t169 * t239) * t221 + (t158 * t231 + t160 * t230) * t197 + (t157 * t231 + t159 * t230) * t196 + (t320 * t223 + t322 * t224) * t179 + (t321 * t223 + t323 * t224) * t178 + (t193 * t231 + t194 * t230 + t210 * t242 + t213 * t239 + t318 * t223 - t319 * t224 + Icges(2,3)) * t226) * t226 / 0.2e1 + V_base(4) * t226 * (Icges(2,5) * t243 - Icges(2,6) * t240) + V_base(5) * t226 * (Icges(2,5) * t240 + Icges(2,6) * t243) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
