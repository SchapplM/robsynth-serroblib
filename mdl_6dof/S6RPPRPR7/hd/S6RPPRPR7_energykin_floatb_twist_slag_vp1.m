% Calculate kinetic energy for
% S6RPPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 01:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRPR7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR7_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPPRPR7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR7_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR7_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRPR7_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:52:47
% EndTime: 2019-03-09 01:52:51
% DurationCPUTime: 3.62s
% Computational Cost: add. (1467->340), mult. (1676->469), div. (0->0), fcn. (1514->10), ass. (0->170)
t309 = Icges(2,4) + Icges(3,6);
t308 = Icges(2,1) + Icges(3,2);
t307 = -Icges(3,4) + Icges(2,5);
t306 = Icges(3,5) - Icges(2,6);
t305 = Icges(2,2) + Icges(3,3);
t228 = cos(qJ(1));
t304 = t309 * t228;
t227 = sin(qJ(1));
t303 = t309 * t227;
t302 = t227 * t308 + t304;
t301 = t228 * t308 - t303;
t300 = t227 * t307 - t228 * t306;
t299 = t227 * t306 + t228 * t307;
t224 = cos(pkin(9));
t222 = sin(pkin(9));
t286 = Icges(4,4) * t222;
t250 = Icges(4,2) * t224 + t286;
t150 = Icges(4,6) * t227 - t228 * t250;
t285 = Icges(4,4) * t224;
t252 = Icges(4,1) * t222 + t285;
t152 = Icges(4,5) * t227 - t228 * t252;
t298 = t150 * t224 + t152 * t222 - t228 * t305 - t303;
t149 = Icges(4,6) * t228 + t227 * t250;
t151 = Icges(4,5) * t228 + t227 * t252;
t297 = t149 * t224 + t151 * t222 + t227 * t305 - t304;
t220 = pkin(9) + qJ(4);
t207 = sin(t220);
t209 = cos(t220);
t223 = cos(pkin(10));
t288 = pkin(5) * t223;
t296 = -pkin(8) * t209 + t207 * t288;
t284 = Icges(5,4) * t207;
t249 = Icges(5,2) * t209 + t284;
t134 = Icges(5,6) * t228 + t227 * t249;
t135 = Icges(5,6) * t227 - t228 * t249;
t283 = Icges(5,4) * t209;
t251 = Icges(5,1) * t207 + t283;
t136 = Icges(5,5) * t228 + t227 * t251;
t137 = Icges(5,5) * t227 - t228 * t251;
t167 = -Icges(5,2) * t207 + t283;
t168 = Icges(5,1) * t209 - t284;
t200 = qJD(4) * t227 + V_base(5);
t201 = qJD(4) * t228 + V_base(4);
t210 = V_base(6) + qJD(1);
t295 = (t134 * t209 + t136 * t207) * t201 + (t135 * t209 + t137 * t207) * t200 + (t167 * t209 + t168 * t207) * t210;
t294 = -pkin(2) - pkin(6);
t290 = pkin(3) * t222;
t289 = pkin(3) * t224;
t280 = qJ(3) * t228;
t219 = pkin(10) + qJ(6);
t206 = sin(t219);
t279 = t206 * t228;
t208 = cos(t219);
t278 = t208 * t228;
t277 = t209 * t227;
t276 = t209 * t228;
t221 = sin(pkin(10));
t275 = t221 * t228;
t274 = t223 * t228;
t273 = t227 * qJ(3);
t272 = t227 * t206;
t271 = t227 * t208;
t270 = t227 * t221;
t269 = t227 * t223;
t266 = qJD(5) * t209;
t265 = qJD(6) * t209;
t193 = t227 * pkin(1) - qJ(2) * t228;
t264 = V_base(4) * t193 + V_base(3);
t263 = V_base(5) * pkin(6) + V_base(1);
t260 = V_base(4) * t273 + t264;
t259 = qJD(2) * t227 + t263;
t258 = -t193 - t273;
t196 = pkin(1) * t228 + t227 * qJ(2);
t257 = -t196 - t280;
t157 = pkin(7) * t227 - t228 * t290;
t256 = -t157 + t258;
t255 = rSges(4,1) * t222 + rSges(4,2) * t224;
t254 = rSges(5,1) * t207 + rSges(5,2) * t209;
t253 = pkin(4) * t207 - qJ(5) * t209;
t248 = Icges(4,5) * t222 + Icges(4,6) * t224;
t247 = Icges(5,5) * t207 + Icges(5,6) * t209;
t177 = -Icges(4,2) * t222 + t285;
t178 = Icges(4,1) * t224 - t286;
t241 = t177 * t224 + t178 * t222;
t240 = -qJD(2) * t228 + t210 * t196 + V_base(2);
t239 = V_base(5) * pkin(2) + qJD(3) * t228 + t259;
t156 = t253 * t228;
t238 = t156 + t256;
t237 = V_base(5) * t289 + t239;
t236 = qJD(3) * t227 + t210 * t280 + t240;
t235 = (Icges(5,3) * t228 + t227 * t247) * t201 + (Icges(5,3) * t227 - t228 * t247) * t200 + (Icges(5,5) * t209 - Icges(5,6) * t207) * t210;
t234 = (Icges(4,3) * t228 + t227 * t248) * V_base(4) + (Icges(4,3) * t227 - t228 * t248) * V_base(5) + (Icges(4,5) * t224 - Icges(4,6) * t222) * t210;
t169 = pkin(4) * t209 + qJ(5) * t207;
t233 = t200 * t169 - t227 * t266 + t237;
t158 = pkin(7) * t228 + t227 * t290;
t232 = V_base(4) * t157 + (-t158 + t257) * V_base(5) + t260;
t231 = qJD(5) * t207 - t201 * t156 + t232;
t230 = t210 * t158 + (-t289 + t294) * V_base(4) + t236;
t155 = t253 * t227;
t229 = t210 * t155 + t228 * t266 + t230;
t198 = rSges(2,1) * t228 - t227 * rSges(2,2);
t197 = -rSges(3,2) * t228 + t227 * rSges(3,3);
t195 = t227 * rSges(2,1) + rSges(2,2) * t228;
t194 = -t227 * rSges(3,2) - rSges(3,3) * t228;
t179 = rSges(4,1) * t224 - rSges(4,2) * t222;
t175 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t174 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t173 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t172 = qJD(6) * t207 + t210;
t170 = rSges(5,1) * t209 - rSges(5,2) * t207;
t164 = -t227 * t265 + t201;
t163 = t228 * t265 + t200;
t162 = -t207 * t274 + t270;
t161 = t207 * t275 + t269;
t160 = t207 * t269 + t275;
t159 = -t207 * t270 + t274;
t154 = t227 * rSges(4,3) - t228 * t255;
t153 = rSges(4,3) * t228 + t227 * t255;
t145 = -t207 * t278 + t272;
t144 = t207 * t279 + t271;
t143 = t207 * t271 + t279;
t142 = -t207 * t272 + t278;
t140 = t227 * rSges(5,3) - t228 * t254;
t139 = rSges(5,3) * t228 + t227 * t254;
t130 = V_base(5) * rSges(2,3) - t195 * t210 + t263;
t129 = t198 * t210 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t128 = rSges(6,3) * t207 + (rSges(6,1) * t223 - rSges(6,2) * t221) * t209;
t127 = Icges(6,5) * t207 + (Icges(6,1) * t223 - Icges(6,4) * t221) * t209;
t126 = Icges(6,6) * t207 + (Icges(6,4) * t223 - Icges(6,2) * t221) * t209;
t125 = Icges(6,3) * t207 + (Icges(6,5) * t223 - Icges(6,6) * t221) * t209;
t123 = t195 * V_base(4) - t198 * V_base(5) + V_base(3);
t122 = rSges(7,3) * t207 + (rSges(7,1) * t208 - rSges(7,2) * t206) * t209;
t121 = Icges(7,5) * t207 + (Icges(7,1) * t208 - Icges(7,4) * t206) * t209;
t120 = Icges(7,6) * t207 + (Icges(7,4) * t208 - Icges(7,2) * t206) * t209;
t119 = Icges(7,3) * t207 + (Icges(7,5) * t208 - Icges(7,6) * t206) * t209;
t118 = pkin(8) * t207 + t209 * t288;
t117 = V_base(5) * rSges(3,1) + (-t193 - t194) * t210 + t259;
t116 = t210 * t197 + (-rSges(3,1) - pkin(6)) * V_base(4) + t240;
t115 = t162 * rSges(6,1) + t161 * rSges(6,2) + rSges(6,3) * t276;
t114 = rSges(6,1) * t160 + rSges(6,2) * t159 - rSges(6,3) * t277;
t113 = Icges(6,1) * t162 + Icges(6,4) * t161 + Icges(6,5) * t276;
t112 = Icges(6,1) * t160 + Icges(6,4) * t159 - Icges(6,5) * t277;
t111 = Icges(6,4) * t162 + Icges(6,2) * t161 + Icges(6,6) * t276;
t110 = Icges(6,4) * t160 + Icges(6,2) * t159 - Icges(6,6) * t277;
t109 = Icges(6,5) * t162 + Icges(6,6) * t161 + Icges(6,3) * t276;
t108 = Icges(6,5) * t160 + Icges(6,6) * t159 - Icges(6,3) * t277;
t107 = pkin(5) * t270 - t296 * t228;
t106 = pkin(5) * t275 + t296 * t227;
t105 = t194 * V_base(4) + (-t196 - t197) * V_base(5) + t264;
t104 = t145 * rSges(7,1) + t144 * rSges(7,2) + rSges(7,3) * t276;
t103 = rSges(7,1) * t143 + rSges(7,2) * t142 - rSges(7,3) * t277;
t102 = Icges(7,1) * t145 + Icges(7,4) * t144 + Icges(7,5) * t276;
t101 = Icges(7,1) * t143 + Icges(7,4) * t142 - Icges(7,5) * t277;
t100 = Icges(7,4) * t145 + Icges(7,2) * t144 + Icges(7,6) * t276;
t99 = Icges(7,4) * t143 + Icges(7,2) * t142 - Icges(7,6) * t277;
t98 = Icges(7,5) * t145 + Icges(7,6) * t144 + Icges(7,3) * t276;
t97 = Icges(7,5) * t143 + Icges(7,6) * t142 - Icges(7,3) * t277;
t96 = t179 * V_base(5) + (-t154 + t258) * t210 + t239;
t95 = t210 * t153 + (-t179 + t294) * V_base(4) + t236;
t94 = V_base(4) * t154 + (-t153 + t257) * V_base(5) + t260;
t93 = t170 * t200 + (-t140 + t256) * t210 + t237;
t92 = t210 * t139 - t201 * t170 + t230;
t91 = -t200 * t139 + t201 * t140 + t232;
t90 = t128 * t200 + (-t115 + t238) * t210 + t233;
t89 = t210 * t114 + (-t128 - t169) * t201 + t229;
t88 = t201 * t115 + (-t114 - t155) * t200 + t231;
t87 = -t104 * t172 + t118 * t200 + t122 * t163 + (-t107 + t238) * t210 + t233;
t86 = t172 * t103 + t210 * t106 - t164 * t122 + (-t118 - t169) * t201 + t229;
t85 = -t163 * t103 + t164 * t104 + t201 * t107 + (-t106 - t155) * t200 + t231;
t1 = t172 * ((t119 * t172 + t163 * t98 + t164 * t97) * t207 + ((t101 * t208 - t206 * t99) * t164 + (-t100 * t206 + t102 * t208) * t163 + (-t120 * t206 + t121 * t208) * t172) * t209) / 0.2e1 + t163 * ((t145 * t101 + t144 * t99 + t97 * t276) * t164 + (t144 * t100 + t145 * t102 + t98 * t276) * t163 + (t119 * t276 + t144 * t120 + t145 * t121) * t172) / 0.2e1 + t164 * ((t143 * t101 + t142 * t99 - t97 * t277) * t164 + (t100 * t142 + t102 * t143 - t98 * t277) * t163 + (-t119 * t277 + t120 * t142 + t121 * t143) * t172) / 0.2e1 + m(1) * (t173 ^ 2 + t174 ^ 2 + t175 ^ 2) / 0.2e1 + m(2) * (t123 ^ 2 + t129 ^ 2 + t130 ^ 2) / 0.2e1 + m(3) * (t105 ^ 2 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + m(5) * (t91 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + m(4) * (t94 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + m(7) * (t85 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + m(6) * (t88 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + ((t108 * t276 + t161 * t110 + t162 * t112) * t201 + (t109 * t276 + t161 * t111 + t162 * t113) * t200 + (t125 * t276 + t161 * t126 + t162 * t127) * t210 + t235 * t227 - t295 * t228) * t200 / 0.2e1 + (t235 * t228 + t295 * t227 + (-t108 * t277 + t159 * t110 + t160 * t112) * t201 + (-t109 * t277 + t111 * t159 + t113 * t160) * t200 + (-t125 * t277 + t126 * t159 + t127 * t160) * t210) * t201 / 0.2e1 + (t234 * t228 + (t241 * t227 + t299) * t210 + (t298 * t227 + t302 * t228 + Icges(1,4)) * V_base(5) + (t297 * t227 + t301 * t228 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t234 * t227 + (-t241 * t228 + t300) * t210 + (t302 * t227 - t298 * t228 + Icges(1,2)) * V_base(5) + (t301 * t227 - t297 * t228 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t108 * t201 + t109 * t200) * t207 + ((-t110 * t221 + t112 * t223) * t201 + (-t111 * t221 + t113 * t223) * t200) * t209 + (-t134 * t207 + t136 * t209) * t201 + (-t135 * t207 + t137 * t209) * t200 + (-t150 * t222 + t152 * t224 + t300) * V_base(5) + (-t149 * t222 + t151 * t224 + t299) * V_base(4) + (-t222 * t177 + t224 * t178 + Icges(3,1) + Icges(2,3) + (-t126 * t221 + t127 * t223 + t168) * t209 + (t125 - t167) * t207) * t210) * t210 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
