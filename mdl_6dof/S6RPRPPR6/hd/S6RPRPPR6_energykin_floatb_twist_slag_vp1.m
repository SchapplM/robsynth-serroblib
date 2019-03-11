% Calculate kinetic energy for
% S6RPRPPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
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
% Datum: 2019-03-09 02:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPPR6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR6_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRPPR6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR6_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR6_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPPR6_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:53:11
% EndTime: 2019-03-09 02:53:14
% DurationCPUTime: 3.10s
% Computational Cost: add. (1489->334), mult. (1698->452), div. (0->0), fcn. (1536->10), ass. (0->166)
t315 = Icges(2,4) + Icges(3,6);
t314 = Icges(2,1) + Icges(3,2);
t313 = -Icges(3,4) + Icges(2,5);
t312 = Icges(3,5) - Icges(2,6);
t311 = Icges(2,2) + Icges(3,3);
t310 = Icges(4,3) + Icges(5,3);
t220 = qJ(3) + pkin(9);
t207 = sin(t220);
t209 = cos(t220);
t225 = sin(qJ(3));
t227 = cos(qJ(3));
t309 = Icges(4,5) * t225 + Icges(5,5) * t207 + Icges(4,6) * t227 + Icges(5,6) * t209;
t228 = cos(qJ(1));
t308 = t315 * t228;
t226 = sin(qJ(1));
t307 = t315 * t226;
t282 = Icges(5,4) * t207;
t249 = Icges(5,2) * t209 + t282;
t135 = Icges(5,6) * t228 + t226 * t249;
t136 = Icges(5,6) * t226 - t228 * t249;
t281 = Icges(5,4) * t209;
t251 = Icges(5,1) * t207 + t281;
t137 = Icges(5,5) * t228 + t226 * t251;
t138 = Icges(5,5) * t226 - t228 * t251;
t284 = Icges(4,4) * t225;
t250 = Icges(4,2) * t227 + t284;
t149 = Icges(4,6) * t228 + t226 * t250;
t150 = Icges(4,6) * t226 - t228 * t250;
t283 = Icges(4,4) * t227;
t252 = Icges(4,1) * t225 + t283;
t151 = Icges(4,5) * t228 + t226 * t252;
t152 = Icges(4,5) * t226 - t228 * t252;
t167 = -Icges(5,2) * t207 + t281;
t168 = Icges(5,1) * t209 - t282;
t184 = -Icges(4,2) * t225 + t283;
t189 = Icges(4,1) * t227 - t284;
t201 = qJD(3) * t226 + V_base(5);
t202 = qJD(3) * t228 + V_base(4);
t210 = V_base(6) + qJD(1);
t306 = t201 * (t136 * t209 + t138 * t207 + t150 * t227 + t152 * t225) + t202 * (t135 * t209 + t137 * t207 + t149 * t227 + t151 * t225) + t210 * (t167 * t209 + t168 * t207 + t184 * t227 + t189 * t225);
t305 = -t311 * t228 - t307;
t304 = t311 * t226 - t308;
t303 = t314 * t226 + t308;
t302 = t314 * t228 - t307;
t299 = (Icges(4,5) * t227 + Icges(5,5) * t209 - Icges(4,6) * t225 - Icges(5,6) * t207) * t210 + (t309 * t226 + t310 * t228) * t202 + (t310 * t226 - t309 * t228) * t201;
t222 = cos(pkin(10));
t287 = pkin(5) * t222;
t295 = -pkin(8) * t209 + t207 * t287;
t291 = pkin(3) * t225;
t290 = pkin(3) * t227;
t289 = pkin(7) * t228;
t288 = t226 * pkin(7);
t219 = pkin(10) + qJ(6);
t206 = sin(t219);
t278 = t206 * t228;
t208 = cos(t219);
t277 = t208 * t228;
t276 = t209 * t226;
t275 = t209 * t228;
t221 = sin(pkin(10));
t274 = t221 * t228;
t273 = t222 * t228;
t272 = t226 * t206;
t271 = t226 * t208;
t270 = t226 * t221;
t269 = t226 * t222;
t253 = pkin(4) * t207 - qJ(5) * t209;
t153 = t253 * t226;
t162 = qJ(4) * t228 + t226 * t291;
t267 = -t153 - t162;
t266 = qJD(5) * t209;
t265 = qJD(6) * t209;
t193 = t226 * pkin(1) - qJ(2) * t228;
t264 = V_base(4) * t193 + V_base(3);
t263 = V_base(5) * pkin(6) + V_base(1);
t169 = pkin(4) * t209 + qJ(5) * t207;
t260 = -t169 - t290;
t259 = -t193 - t288;
t258 = qJD(2) * t226 + t263;
t161 = qJ(4) * t226 - t228 * t291;
t257 = -t161 + t259;
t256 = V_base(5) * pkin(2) + t258;
t255 = rSges(4,1) * t225 + rSges(4,2) * t227;
t254 = rSges(5,1) * t207 + rSges(5,2) * t209;
t197 = pkin(1) * t228 + t226 * qJ(2);
t240 = -qJD(2) * t228 + t210 * t197 + V_base(2);
t154 = t253 * t228;
t239 = t154 + t257;
t238 = qJD(4) * t228 + t201 * t290 + t256;
t235 = V_base(4) * t288 + (-t197 - t289) * V_base(5) + t264;
t234 = t210 * t289 + (-pkin(2) - pkin(6)) * V_base(4) + t240;
t233 = t202 * t161 + t235;
t232 = t201 * t169 - t226 * t266 + t238;
t231 = qJD(4) * t226 + t210 * t162 + t234;
t230 = qJD(5) * t207 - t202 * t154 + t233;
t229 = t210 * t153 + t228 * t266 + t231;
t199 = rSges(2,1) * t228 - t226 * rSges(2,2);
t198 = -rSges(3,2) * t228 + t226 * rSges(3,3);
t196 = rSges(4,1) * t227 - rSges(4,2) * t225;
t195 = t226 * rSges(2,1) + rSges(2,2) * t228;
t194 = -t226 * rSges(3,2) - rSges(3,3) * t228;
t176 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t175 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t174 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t173 = qJD(6) * t207 + t210;
t170 = rSges(5,1) * t209 - rSges(5,2) * t207;
t164 = -t226 * t265 + t202;
t163 = t228 * t265 + t201;
t160 = -t207 * t273 + t270;
t159 = t207 * t274 + t269;
t158 = t207 * t269 + t274;
t157 = -t207 * t270 + t273;
t156 = t226 * rSges(4,3) - t228 * t255;
t155 = rSges(4,3) * t228 + t226 * t255;
t146 = -t207 * t277 + t272;
t145 = t207 * t278 + t271;
t144 = t207 * t271 + t278;
t143 = -t207 * t272 + t277;
t140 = t226 * rSges(5,3) - t228 * t254;
t139 = rSges(5,3) * t228 + t226 * t254;
t131 = V_base(5) * rSges(2,3) - t195 * t210 + t263;
t130 = t199 * t210 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t128 = rSges(6,3) * t207 + (rSges(6,1) * t222 - rSges(6,2) * t221) * t209;
t127 = Icges(6,5) * t207 + (Icges(6,1) * t222 - Icges(6,4) * t221) * t209;
t126 = Icges(6,6) * t207 + (Icges(6,4) * t222 - Icges(6,2) * t221) * t209;
t125 = Icges(6,3) * t207 + (Icges(6,5) * t222 - Icges(6,6) * t221) * t209;
t123 = t195 * V_base(4) - t199 * V_base(5) + V_base(3);
t122 = rSges(7,3) * t207 + (rSges(7,1) * t208 - rSges(7,2) * t206) * t209;
t121 = Icges(7,5) * t207 + (Icges(7,1) * t208 - Icges(7,4) * t206) * t209;
t120 = Icges(7,6) * t207 + (Icges(7,4) * t208 - Icges(7,2) * t206) * t209;
t119 = Icges(7,3) * t207 + (Icges(7,5) * t208 - Icges(7,6) * t206) * t209;
t118 = pkin(8) * t207 + t209 * t287;
t117 = V_base(5) * rSges(3,1) + (-t193 - t194) * t210 + t258;
t116 = t210 * t198 + (-rSges(3,1) - pkin(6)) * V_base(4) + t240;
t115 = t160 * rSges(6,1) + t159 * rSges(6,2) + rSges(6,3) * t275;
t114 = rSges(6,1) * t158 + rSges(6,2) * t157 - rSges(6,3) * t276;
t113 = Icges(6,1) * t160 + Icges(6,4) * t159 + Icges(6,5) * t275;
t112 = Icges(6,1) * t158 + Icges(6,4) * t157 - Icges(6,5) * t276;
t111 = Icges(6,4) * t160 + Icges(6,2) * t159 + Icges(6,6) * t275;
t110 = Icges(6,4) * t158 + Icges(6,2) * t157 - Icges(6,6) * t276;
t109 = Icges(6,5) * t160 + Icges(6,6) * t159 + Icges(6,3) * t275;
t108 = Icges(6,5) * t158 + Icges(6,6) * t157 - Icges(6,3) * t276;
t107 = pkin(5) * t270 - t228 * t295;
t106 = pkin(5) * t274 + t226 * t295;
t105 = t194 * V_base(4) + (-t197 - t198) * V_base(5) + t264;
t104 = t146 * rSges(7,1) + t145 * rSges(7,2) + rSges(7,3) * t275;
t103 = rSges(7,1) * t144 + rSges(7,2) * t143 - rSges(7,3) * t276;
t102 = Icges(7,1) * t146 + Icges(7,4) * t145 + Icges(7,5) * t275;
t101 = Icges(7,1) * t144 + Icges(7,4) * t143 - Icges(7,5) * t276;
t100 = Icges(7,4) * t146 + Icges(7,2) * t145 + Icges(7,6) * t275;
t99 = Icges(7,4) * t144 + Icges(7,2) * t143 - Icges(7,6) * t276;
t98 = Icges(7,5) * t146 + Icges(7,6) * t145 + Icges(7,3) * t275;
t97 = Icges(7,5) * t144 + Icges(7,6) * t143 - Icges(7,3) * t276;
t96 = t196 * t201 + (-t156 + t259) * t210 + t256;
t95 = t210 * t155 - t202 * t196 + t234;
t94 = -t201 * t155 + t202 * t156 + t235;
t93 = t170 * t201 + (-t140 + t257) * t210 + t238;
t92 = t210 * t139 + (-t170 - t290) * t202 + t231;
t91 = t202 * t140 + (-t139 - t162) * t201 + t233;
t90 = t128 * t201 + (-t115 + t239) * t210 + t232;
t89 = t210 * t114 + (-t128 + t260) * t202 + t229;
t88 = t202 * t115 + (-t114 + t267) * t201 + t230;
t87 = -t104 * t173 + t118 * t201 + t122 * t163 + (-t107 + t239) * t210 + t232;
t86 = t173 * t103 + t210 * t106 - t164 * t122 + (-t118 + t260) * t202 + t229;
t85 = -t163 * t103 + t164 * t104 + t202 * t107 + (-t106 + t267) * t201 + t230;
t1 = t164 * ((t101 * t144 + t143 * t99 - t97 * t276) * t164 + (t100 * t143 + t102 * t144 - t276 * t98) * t163 + (-t119 * t276 + t120 * t143 + t121 * t144) * t173) / 0.2e1 + t163 * ((t146 * t101 + t145 * t99 + t275 * t97) * t164 + (t145 * t100 + t146 * t102 + t275 * t98) * t163 + (t119 * t275 + t145 * t120 + t146 * t121) * t173) / 0.2e1 + t173 * ((t119 * t173 + t98 * t163 + t97 * t164) * t207 + ((t101 * t208 - t206 * t99) * t164 + (-t100 * t206 + t102 * t208) * t163 + (-t120 * t206 + t121 * t208) * t173) * t209) / 0.2e1 + m(7) * (t85 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + m(6) * (t88 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + m(5) * (t91 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + m(4) * (t94 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + m(3) * (t105 ^ 2 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + m(2) * (t123 ^ 2 + t130 ^ 2 + t131 ^ 2) / 0.2e1 + m(1) * (t174 ^ 2 + t175 ^ 2 + t176 ^ 2) / 0.2e1 + ((t108 * t275 + t159 * t110 + t160 * t112) * t202 + (t109 * t275 + t159 * t111 + t160 * t113) * t201 + (t125 * t275 + t159 * t126 + t160 * t127) * t210 - t306 * t228 + t299 * t226) * t201 / 0.2e1 + ((-t108 * t276 + t110 * t157 + t112 * t158) * t202 + (-t109 * t276 + t111 * t157 + t113 * t158) * t201 + (-t125 * t276 + t126 * t157 + t127 * t158) * t210 + t299 * t228 + t306 * t226) * t202 / 0.2e1 + ((t226 * t305 + t303 * t228 + Icges(1,4)) * V_base(5) + (t304 * t226 + t302 * t228 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t303 * t226 - t305 * t228 + Icges(1,2)) * V_base(5) + (t226 * t302 - t228 * t304 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t149 * t225 + t151 * t227 + (-t110 * t221 + t112 * t222 + t137) * t209 + (t108 - t135) * t207) * t202 + (-t150 * t225 + t152 * t227 + (-t111 * t221 + t113 * t222 + t138) * t209 + (t109 - t136) * t207) * t201 + (-t184 * t225 + t189 * t227 + Icges(3,1) + Icges(2,3) + (-t126 * t221 + t127 * t222 + t168) * t209 + (t125 - t167) * t207) * t210) * t210 / 0.2e1 + t210 * V_base(5) * (t313 * t226 - t312 * t228) + t210 * V_base(4) * (t312 * t226 + t313 * t228) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
