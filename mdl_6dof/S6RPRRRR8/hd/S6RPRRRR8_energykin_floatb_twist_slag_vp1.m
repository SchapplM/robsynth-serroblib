% Calculate kinetic energy for
% S6RPRRRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% Datum: 2019-03-09 07:22
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRR8_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR8_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR8_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRRR8_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR8_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR8_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR8_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRR8_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:19:41
% EndTime: 2019-03-09 07:19:44
% DurationCPUTime: 3.28s
% Computational Cost: add. (1603->334), mult. (1800->490), div. (0->0), fcn. (1638->10), ass. (0->169)
t309 = Icges(2,4) + Icges(3,6);
t308 = Icges(2,1) + Icges(3,2);
t307 = -Icges(3,4) + Icges(2,5);
t306 = Icges(3,5) - Icges(2,6);
t305 = Icges(2,2) + Icges(3,3);
t231 = cos(qJ(1));
t304 = t309 * t231;
t228 = sin(qJ(1));
t303 = t309 * t228;
t302 = -t231 * t305 - t303;
t301 = t228 * t305 - t304;
t300 = t228 * t308 + t304;
t299 = t231 * t308 - t303;
t225 = qJ(3) + qJ(4);
t218 = sin(t225);
t220 = cos(t225);
t229 = cos(qJ(5));
t288 = pkin(5) * t229;
t296 = -pkin(10) * t220 + t218 * t288;
t282 = Icges(5,4) * t218;
t253 = Icges(5,2) * t220 + t282;
t138 = Icges(5,6) * t231 + t228 * t253;
t139 = Icges(5,6) * t228 - t231 * t253;
t281 = Icges(5,4) * t220;
t255 = Icges(5,1) * t218 + t281;
t140 = Icges(5,5) * t231 + t228 * t255;
t141 = Icges(5,5) * t228 - t231 * t255;
t170 = -Icges(5,2) * t218 + t281;
t171 = Icges(5,1) * t220 - t282;
t206 = qJD(3) * t228 + V_base(5);
t174 = qJD(4) * t228 + t206;
t207 = qJD(3) * t231 + V_base(4);
t175 = qJD(4) * t231 + t207;
t211 = V_base(6) + qJD(1);
t295 = (t138 * t220 + t140 * t218) * t175 + (t139 * t220 + t141 * t218) * t174 + (t170 * t220 + t171 * t218) * t211;
t230 = cos(qJ(3));
t227 = sin(qJ(3));
t284 = Icges(4,4) * t227;
t254 = Icges(4,2) * t230 + t284;
t153 = Icges(4,6) * t231 + t228 * t254;
t154 = Icges(4,6) * t228 - t231 * t254;
t283 = Icges(4,4) * t230;
t256 = Icges(4,1) * t227 + t283;
t155 = Icges(4,5) * t231 + t228 * t256;
t156 = Icges(4,5) * t228 - t231 * t256;
t189 = -Icges(4,2) * t227 + t283;
t194 = Icges(4,1) * t230 - t284;
t294 = (t153 * t230 + t155 * t227) * t207 + (t154 * t230 + t156 * t227) * t206 + (t189 * t230 + t194 * t227) * t211;
t292 = pkin(3) * t227;
t291 = pkin(3) * t230;
t290 = t228 * pkin(7);
t289 = t231 * pkin(7);
t224 = qJ(5) + qJ(6);
t217 = sin(t224);
t278 = t217 * t228;
t277 = t217 * t231;
t219 = cos(t224);
t276 = t219 * t228;
t275 = t219 * t231;
t274 = t220 * t228;
t273 = t220 * t231;
t226 = sin(qJ(5));
t272 = t226 * t228;
t271 = t226 * t231;
t270 = t228 * t229;
t269 = t229 * t231;
t268 = qJD(5) * t220;
t197 = pkin(1) * t228 - qJ(2) * t231;
t267 = V_base(4) * t197 + V_base(3);
t266 = V_base(5) * pkin(6) + V_base(1);
t263 = -t197 - t290;
t262 = qJD(2) * t228 + t266;
t181 = qJD(5) * t218 + t211;
t162 = pkin(8) * t228 - t231 * t292;
t261 = -t162 + t263;
t260 = V_base(5) * pkin(2) + t262;
t259 = pkin(4) * t218 - pkin(9) * t220;
t258 = rSges(4,1) * t227 + rSges(4,2) * t230;
t257 = rSges(5,1) * t218 + rSges(5,2) * t220;
t145 = t231 * t268 + t174;
t252 = Icges(4,5) * t227 + Icges(4,6) * t230;
t251 = Icges(5,5) * t218 + Icges(5,6) * t220;
t201 = pkin(1) * t231 + qJ(2) * t228;
t244 = -qJD(2) * t231 + t211 * t201 + V_base(2);
t243 = t206 * t291 + t260;
t242 = (Icges(5,3) * t231 + t228 * t251) * t175 + (Icges(5,3) * t228 - t231 * t251) * t174 + (Icges(5,5) * t220 - Icges(5,6) * t218) * t211;
t241 = (Icges(4,3) * t231 + t228 * t252) * t207 + (Icges(4,3) * t228 - t231 * t252) * t206 + (Icges(4,5) * t230 - Icges(4,6) * t227) * t211;
t240 = V_base(4) * t290 + (-t201 - t289) * V_base(5) + t267;
t239 = t211 * t289 + (-pkin(2) - pkin(6)) * V_base(4) + t244;
t163 = pkin(8) * t231 + t228 * t292;
t238 = t207 * t162 - t163 * t206 + t240;
t160 = t259 * t231;
t173 = pkin(4) * t220 + pkin(9) * t218;
t237 = t174 * t173 + (t160 + t261) * t211 + t243;
t236 = t211 * t163 - t207 * t291 + t239;
t159 = t259 * t228;
t235 = -t159 * t174 - t175 * t160 + t238;
t234 = t211 * t159 - t173 * t175 + t236;
t203 = rSges(2,1) * t231 - rSges(2,2) * t228;
t202 = -rSges(3,2) * t231 + rSges(3,3) * t228;
t200 = rSges(4,1) * t230 - rSges(4,2) * t227;
t199 = rSges(2,1) * t228 + rSges(2,2) * t231;
t198 = -rSges(3,2) * t228 - rSges(3,3) * t231;
t180 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t179 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t178 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t172 = rSges(5,1) * t220 - rSges(5,2) * t218;
t167 = -t218 * t269 + t272;
t166 = t218 * t271 + t270;
t165 = t218 * t270 + t271;
t164 = -t218 * t272 + t269;
t161 = qJD(6) * t218 + t181;
t158 = rSges(4,3) * t228 - t231 * t258;
t157 = rSges(4,3) * t231 + t228 * t258;
t150 = -t218 * t275 + t278;
t149 = t218 * t277 + t276;
t148 = t218 * t276 + t277;
t147 = -t218 * t278 + t275;
t146 = -t228 * t268 + t175;
t143 = rSges(5,3) * t228 - t231 * t257;
t142 = rSges(5,3) * t231 + t228 * t257;
t134 = rSges(6,3) * t218 + (rSges(6,1) * t229 - rSges(6,2) * t226) * t220;
t133 = Icges(6,5) * t218 + (Icges(6,1) * t229 - Icges(6,4) * t226) * t220;
t132 = Icges(6,6) * t218 + (Icges(6,4) * t229 - Icges(6,2) * t226) * t220;
t131 = Icges(6,3) * t218 + (Icges(6,5) * t229 - Icges(6,6) * t226) * t220;
t130 = V_base(5) * rSges(2,3) - t199 * t211 + t266;
t129 = t203 * t211 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t126 = t199 * V_base(4) - t203 * V_base(5) + V_base(3);
t125 = rSges(7,3) * t218 + (rSges(7,1) * t219 - rSges(7,2) * t217) * t220;
t124 = Icges(7,5) * t218 + (Icges(7,1) * t219 - Icges(7,4) * t217) * t220;
t123 = Icges(7,6) * t218 + (Icges(7,4) * t219 - Icges(7,2) * t217) * t220;
t122 = Icges(7,3) * t218 + (Icges(7,5) * t219 - Icges(7,6) * t217) * t220;
t121 = (-qJD(5) - qJD(6)) * t274 + t175;
t120 = qJD(6) * t273 + t145;
t118 = pkin(10) * t218 + t220 * t288;
t117 = V_base(5) * rSges(3,1) + (-t197 - t198) * t211 + t262;
t116 = t202 * t211 + (-rSges(3,1) - pkin(6)) * V_base(4) + t244;
t115 = rSges(6,1) * t167 + rSges(6,2) * t166 + rSges(6,3) * t273;
t114 = rSges(6,1) * t165 + rSges(6,2) * t164 - rSges(6,3) * t274;
t113 = Icges(6,1) * t167 + Icges(6,4) * t166 + Icges(6,5) * t273;
t112 = Icges(6,1) * t165 + Icges(6,4) * t164 - Icges(6,5) * t274;
t111 = Icges(6,4) * t167 + Icges(6,2) * t166 + Icges(6,6) * t273;
t110 = Icges(6,4) * t165 + Icges(6,2) * t164 - Icges(6,6) * t274;
t109 = Icges(6,5) * t167 + Icges(6,6) * t166 + Icges(6,3) * t273;
t108 = Icges(6,5) * t165 + Icges(6,6) * t164 - Icges(6,3) * t274;
t107 = pkin(5) * t272 - t231 * t296;
t106 = pkin(5) * t271 + t228 * t296;
t105 = t198 * V_base(4) + (-t201 - t202) * V_base(5) + t267;
t104 = rSges(7,1) * t150 + rSges(7,2) * t149 + rSges(7,3) * t273;
t103 = rSges(7,1) * t148 + rSges(7,2) * t147 - rSges(7,3) * t274;
t102 = Icges(7,1) * t150 + Icges(7,4) * t149 + Icges(7,5) * t273;
t101 = Icges(7,1) * t148 + Icges(7,4) * t147 - Icges(7,5) * t274;
t100 = Icges(7,4) * t150 + Icges(7,2) * t149 + Icges(7,6) * t273;
t99 = Icges(7,4) * t148 + Icges(7,2) * t147 - Icges(7,6) * t274;
t98 = Icges(7,5) * t150 + Icges(7,6) * t149 + Icges(7,3) * t273;
t97 = Icges(7,5) * t148 + Icges(7,6) * t147 - Icges(7,3) * t274;
t96 = t200 * t206 + (-t158 + t263) * t211 + t260;
t95 = t157 * t211 - t200 * t207 + t239;
t94 = -t157 * t206 + t158 * t207 + t240;
t93 = t172 * t174 + (-t143 + t261) * t211 + t243;
t92 = t142 * t211 - t172 * t175 + t236;
t91 = -t142 * t174 + t143 * t175 + t238;
t90 = -t115 * t181 + t134 * t145 + t237;
t89 = t114 * t181 - t134 * t146 + t234;
t88 = -t114 * t145 + t115 * t146 + t235;
t87 = -t104 * t161 - t107 * t181 + t118 * t145 + t120 * t125 + t237;
t86 = t103 * t161 + t106 * t181 - t118 * t146 - t121 * t125 + t234;
t85 = -t103 * t120 + t104 * t121 - t106 * t145 + t107 * t146 + t235;
t1 = t145 * ((t108 * t273 + t110 * t166 + t112 * t167) * t146 + (t109 * t273 + t111 * t166 + t113 * t167) * t145 + (t131 * t273 + t132 * t166 + t133 * t167) * t181) / 0.2e1 + t120 * ((t101 * t150 + t149 * t99 + t273 * t97) * t121 + (t100 * t149 + t102 * t150 + t98 * t273) * t120 + (t122 * t273 + t123 * t149 + t124 * t150) * t161) / 0.2e1 + t146 * ((-t108 * t274 + t110 * t164 + t112 * t165) * t146 + (-t109 * t274 + t111 * t164 + t113 * t165) * t145 + (-t131 * t274 + t132 * t164 + t133 * t165) * t181) / 0.2e1 + t121 * ((t101 * t148 + t147 * t99 - t97 * t274) * t121 + (t100 * t147 + t102 * t148 - t274 * t98) * t120 + (-t122 * t274 + t123 * t147 + t124 * t148) * t161) / 0.2e1 + t207 * (t294 * t228 + t241 * t231) / 0.2e1 + t206 * (t241 * t228 - t294 * t231) / 0.2e1 + t175 * (t295 * t228 + t242 * t231) / 0.2e1 + t174 * (t242 * t228 - t295 * t231) / 0.2e1 + m(1) * (t178 ^ 2 + t179 ^ 2 + t180 ^ 2) / 0.2e1 + m(2) * (t126 ^ 2 + t129 ^ 2 + t130 ^ 2) / 0.2e1 + m(3) * (t105 ^ 2 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + m(5) * (t91 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + m(4) * (t94 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + m(7) * (t85 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + m(6) * (t88 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + t161 * ((t98 * t120 + t97 * t121 + t122 * t161) * t218 + ((t101 * t219 - t217 * t99) * t121 + (-t100 * t217 + t102 * t219) * t120 + (-t123 * t217 + t124 * t219) * t161) * t220) / 0.2e1 + t181 * ((t108 * t146 + t109 * t145 + t131 * t181) * t218 + ((-t110 * t226 + t112 * t229) * t146 + (-t111 * t226 + t113 * t229) * t145 + (-t132 * t226 + t133 * t229) * t181) * t220) / 0.2e1 + ((t228 * t302 + t300 * t231 + Icges(1,4)) * V_base(5) + (t301 * t228 + t299 * t231 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t300 * t228 - t302 * t231 + Icges(1,2)) * V_base(5) + (t228 * t299 - t231 * t301 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t153 * t227 + t155 * t230) * t207 + (-t154 * t227 + t156 * t230) * t206 + (-t138 * t218 + t140 * t220) * t175 + (-t139 * t218 + t141 * t220) * t174 + (-t170 * t218 + t171 * t220 - t189 * t227 + t194 * t230 + Icges(3,1) + Icges(2,3)) * t211) * t211 / 0.2e1 + t211 * V_base(5) * (t228 * t307 - t231 * t306) + t211 * V_base(4) * (t228 * t306 + t307 * t231) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
