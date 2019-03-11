% Calculate joint inertia matrix for
% S6RPRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP8_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP8_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP8_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRP8_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:22:57
% EndTime: 2019-03-09 06:23:05
% DurationCPUTime: 3.50s
% Computational Cost: add. (6071->375), mult. (8175->540), div. (0->0), fcn. (8774->8), ass. (0->184)
t189 = qJ(3) + qJ(4);
t180 = sin(t189);
t193 = cos(qJ(5));
t195 = cos(qJ(1));
t238 = t195 * t193;
t190 = sin(qJ(5));
t192 = sin(qJ(1));
t242 = t192 * t190;
t149 = t180 * t242 - t238;
t239 = t195 * t190;
t241 = t192 * t193;
t150 = t180 * t241 + t239;
t181 = cos(t189);
t246 = t181 * t192;
t79 = Icges(7,5) * t150 - Icges(7,6) * t246 + Icges(7,3) * t149;
t85 = Icges(6,4) * t150 - Icges(6,2) * t149 - Icges(6,6) * t246;
t290 = t79 - t85;
t151 = t180 * t239 + t241;
t153 = -t180 * t238 + t242;
t244 = t181 * t195;
t80 = Icges(7,5) * t153 + Icges(7,6) * t244 - Icges(7,3) * t151;
t86 = Icges(6,4) * t153 + Icges(6,2) * t151 + Icges(6,6) * t244;
t289 = t80 - t86;
t81 = Icges(6,5) * t150 - Icges(6,6) * t149 - Icges(6,3) * t246;
t83 = Icges(7,4) * t150 - Icges(7,2) * t246 + Icges(7,6) * t149;
t288 = t81 + t83;
t82 = Icges(6,5) * t153 + Icges(6,6) * t151 + Icges(6,3) * t244;
t84 = Icges(7,4) * t153 + Icges(7,2) * t244 - Icges(7,6) * t151;
t287 = t82 + t84;
t87 = Icges(7,1) * t150 - Icges(7,4) * t246 + Icges(7,5) * t149;
t89 = Icges(6,1) * t150 - Icges(6,4) * t149 - Icges(6,5) * t246;
t286 = t87 + t89;
t88 = Icges(7,1) * t153 + Icges(7,4) * t244 - Icges(7,5) * t151;
t90 = Icges(6,1) * t153 + Icges(6,4) * t151 + Icges(6,5) * t244;
t285 = t88 + t90;
t276 = rSges(7,3) + qJ(6);
t277 = rSges(7,1) + pkin(5);
t284 = -t276 * t151 + t277 * t153;
t283 = t149 * t290 + t286 * t150 - t288 * t246;
t282 = t149 * t289 + t150 * t285 - t246 * t287;
t281 = -t151 * t290 + t286 * t153 + t288 * t244;
t280 = -t151 * t289 + t153 * t285 + t244 * t287;
t113 = Icges(7,6) * t180 + (Icges(7,5) * t193 + Icges(7,3) * t190) * t181;
t115 = Icges(7,2) * t180 + (Icges(7,4) * t193 + Icges(7,6) * t190) * t181;
t117 = Icges(7,4) * t180 + (Icges(7,1) * t193 + Icges(7,5) * t190) * t181;
t54 = t113 * t149 - t115 * t246 + t117 * t150;
t114 = Icges(6,3) * t180 + (Icges(6,5) * t193 - Icges(6,6) * t190) * t181;
t116 = Icges(6,6) * t180 + (Icges(6,4) * t193 - Icges(6,2) * t190) * t181;
t118 = Icges(6,5) * t180 + (Icges(6,1) * t193 - Icges(6,4) * t190) * t181;
t55 = -t114 * t246 - t116 * t149 + t118 * t150;
t279 = t54 + t55;
t56 = -t113 * t151 + t115 * t244 + t117 * t153;
t57 = t114 * t244 + t116 * t151 + t118 * t153;
t278 = t56 + t57;
t275 = (-t192 * t283 + t195 * t282) * t181 + t279 * t180;
t274 = (-t192 * t281 + t195 * t280) * t181 + t278 * t180;
t273 = t192 * t282 + t195 * t283;
t272 = t192 * t280 + t195 * t281;
t271 = rSges(7,2) * t244 + t284;
t236 = t180 * rSges(7,2) + (t276 * t190 + t193 * t277) * t181;
t247 = t181 * t190;
t270 = t113 * t247 + (t117 + t118) * t181 * t193 + (t114 + t115) * t180;
t269 = t192 * t195;
t268 = (rSges(5,1) * t180 + rSges(5,2) * t181) * t195;
t191 = sin(qJ(3));
t194 = cos(qJ(3));
t267 = (rSges(4,1) * t191 + rSges(4,2) * t194) * t195;
t266 = t276 * t149 + t277 * t150;
t187 = t192 ^ 2;
t188 = t195 ^ 2;
t264 = t192 / 0.2e1;
t263 = t195 / 0.2e1;
t262 = pkin(3) * t194;
t261 = (-t116 * t247 + t270) * t180;
t41 = t180 * t83 + (t190 * t79 + t193 * t87) * t181;
t260 = t41 * t195;
t42 = t180 * t84 + (t190 * t80 + t193 * t88) * t181;
t259 = t42 * t192;
t43 = t180 * t81 + (-t190 * t85 + t193 * t89) * t181;
t258 = t43 * t195;
t44 = t180 * t82 + (-t190 * t86 + t193 * t90) * t181;
t257 = t44 * t192;
t256 = -rSges(7,2) * t246 + t266;
t171 = t195 * t180 * pkin(4);
t131 = t195 * (pkin(9) * t244 - t171);
t217 = -t153 * rSges(6,1) - t151 * rSges(6,2);
t94 = rSges(6,3) * t244 - t217;
t254 = t195 * t94 + t131;
t248 = t180 * t192;
t170 = pkin(4) * t248;
t147 = -pkin(9) * t246 + t170;
t234 = t150 * rSges(6,1) - t149 * rSges(6,2);
t92 = -rSges(6,3) * t246 + t234;
t253 = -t147 - t92;
t252 = Icges(4,4) * t191;
t251 = Icges(4,4) * t194;
t250 = Icges(5,4) * t180;
t249 = Icges(5,4) * t181;
t243 = t191 * t192;
t120 = t180 * rSges(6,3) + (rSges(6,1) * t193 - rSges(6,2) * t190) * t181;
t111 = t192 * t120;
t240 = t192 * t194;
t160 = t181 * pkin(4) + t180 * pkin(9);
t154 = t192 * t160;
t96 = t154 + t111;
t196 = -pkin(8) - pkin(7);
t233 = t195 * t191 * pkin(3) + t192 * t196;
t232 = t195 * pkin(1) + t192 * qJ(2);
t231 = t187 + t188;
t230 = t271 * t195 + t131;
t229 = -t147 - t256;
t68 = t236 * t192 + t154;
t129 = rSges(5,1) * t248 + rSges(5,2) * t246 + t195 * rSges(5,3);
t227 = rSges(4,1) * t243 + rSges(4,2) * t240 + t195 * rSges(4,3);
t183 = t195 * qJ(2);
t226 = t183 + t233;
t225 = t181 * (-rSges(7,2) - pkin(9));
t224 = t181 * (-rSges(6,3) - pkin(9));
t221 = -t160 - t262;
t220 = t181 * t236;
t210 = Icges(5,5) * t180 + Icges(5,6) * t181;
t123 = Icges(5,3) * t195 + t192 * t210;
t124 = Icges(5,3) * t192 - t195 * t210;
t216 = (t187 * t124 + t272) * t192 + (t188 * t123 + (t192 * t123 + t195 * t124) * t192 + t273) * t195;
t215 = Icges(4,1) * t191 + t251;
t214 = Icges(5,1) * t180 + t249;
t213 = Icges(4,2) * t194 + t252;
t212 = Icges(5,2) * t181 + t250;
t211 = Icges(4,5) * t191 + Icges(4,6) * t194;
t159 = t181 * rSges(5,1) - t180 * rSges(5,2);
t176 = pkin(3) * t240;
t121 = t159 * t192 + t176;
t122 = (-t159 - t262) * t195;
t209 = t121 * t192 - t122 * t195;
t157 = -Icges(5,2) * t180 + t249;
t158 = Icges(5,1) * t181 - t250;
t204 = t157 * t181 + t158 * t180;
t175 = pkin(3) * t243;
t203 = -t195 * t196 + t175 + t232;
t98 = t268 + (-rSges(5,3) - pkin(1)) * t192 + t226;
t99 = t203 + t129;
t202 = m(5) * (t192 * t98 - t195 * t99);
t201 = -t192 * pkin(1) + t171 + t226;
t106 = t183 + t267 + (-rSges(4,3) - pkin(1) - pkin(7)) * t192;
t107 = t195 * pkin(7) + t227 + t232;
t200 = m(4) * (t106 * t192 - t107 * t195);
t199 = t170 + t203;
t198 = (t259 + t260 + t257 + t258) * t180 / 0.2e1 + t274 * t264 + t275 * t263 - t273 * t246 / 0.2e1 + t272 * t244 / 0.2e1;
t156 = Icges(5,5) * t181 - Icges(5,6) * t180;
t197 = t260 / 0.2e1 + t259 / 0.2e1 + t258 / 0.2e1 + t257 / 0.2e1 + (-(Icges(5,6) * t192 - t195 * t212) * t180 + (Icges(5,5) * t192 - t195 * t214) * t181 + t192 * t156 - t195 * t204 + t278) * t264 + (-(Icges(5,6) * t195 + t192 * t212) * t180 + (Icges(5,5) * t195 + t192 * t214) * t181 + t195 * t156 + t192 * t204 + t279) * t263;
t167 = t195 * rSges(2,1) - t192 * rSges(2,2);
t166 = t194 * rSges(4,1) - t191 * rSges(4,2);
t165 = -t192 * rSges(2,1) - t195 * rSges(2,2);
t148 = t175 + (-pkin(7) - t196) * t195;
t146 = -t195 * rSges(3,2) + t192 * rSges(3,3) + t232;
t145 = t195 * rSges(3,3) + t183 + (rSges(3,2) - pkin(1)) * t192;
t134 = Icges(4,3) * t192 - t195 * t211;
t133 = Icges(4,3) * t195 + t192 * t211;
t132 = t195 * (-t192 * pkin(7) - t233);
t112 = t195 * (t192 * rSges(5,3) - t268);
t97 = (-t120 - t160) * t195;
t95 = -t192 * t227 + (t192 * rSges(4,3) - t267) * t195;
t76 = (-t120 + t221) * t195;
t75 = t176 + t96;
t74 = -t129 * t192 + t112;
t69 = (-t160 - t236) * t195;
t67 = (t221 - t236) * t195;
t66 = t176 + t68;
t65 = t120 * t244 - t180 * t94;
t64 = t111 * t181 + t180 * t92;
t63 = t192 * t224 + t199 + t234;
t62 = t195 * t224 + t201 + t217;
t61 = t112 + t132 + (-t129 - t148) * t192;
t58 = (-t192 * t94 - t195 * t92) * t181;
t53 = t192 * t225 + t199 + t266;
t52 = t195 * t225 + t201 - t284;
t47 = t192 * t253 + t254;
t46 = -t180 * t271 + t195 * t220;
t45 = t180 * t256 + t192 * t220;
t36 = t132 + (-t148 + t253) * t192 + t254;
t25 = (-t192 * t271 - t195 * t256) * t181;
t24 = t192 * t229 + t230;
t23 = t132 + (-t148 + t229) * t192 + t230;
t1 = [-t180 * t157 - t191 * (-Icges(4,2) * t191 + t251) + t194 * (Icges(4,1) * t194 - t252) + Icges(3,1) + Icges(2,3) + (-t116 * t190 + t158) * t181 + m(6) * (t62 ^ 2 + t63 ^ 2) + m(7) * (t52 ^ 2 + t53 ^ 2) + m(5) * (t98 ^ 2 + t99 ^ 2) + m(4) * (t106 ^ 2 + t107 ^ 2) + m(3) * (t145 ^ 2 + t146 ^ 2) + m(2) * (t165 ^ 2 + t167 ^ 2) + t270; m(6) * (t192 * t62 - t195 * t63) + m(7) * (t192 * t52 - t195 * t53) + t202 + t200 + m(3) * (t145 * t192 - t146 * t195); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * t231; m(6) * (t62 * t75 + t63 * t76) + m(7) * (t52 * t66 + t53 * t67) + m(5) * (t121 * t98 + t122 * t99) + t166 * t200 + t197 + (t187 / 0.2e1 + t188 / 0.2e1) * (Icges(4,5) * t194 - Icges(4,6) * t191) + (-t191 * (Icges(4,6) * t192 - t195 * t213) + t194 * (Icges(4,5) * t192 - t195 * t215)) * t264 + (-t191 * (Icges(4,6) * t195 + t192 * t213) + t194 * (Icges(4,5) * t195 + t192 * t215)) * t263; m(5) * t209 + m(6) * (t192 * t75 - t195 * t76) + m(7) * (t192 * t66 - t195 * t67) + m(4) * t231 * t166; m(7) * (t23 ^ 2 + t66 ^ 2 + t67 ^ 2) + m(6) * (t36 ^ 2 + t75 ^ 2 + t76 ^ 2) + m(5) * (t121 ^ 2 + t122 ^ 2 + t61 ^ 2) + t192 * (t133 * t269 + t187 * t134) + t195 * (t188 * t133 + t134 * t269) + m(4) * (t166 ^ 2 * t231 + t95 ^ 2) + t216; m(6) * (t62 * t96 + t63 * t97) + m(7) * (t52 * t68 + t53 * t69) + t159 * t202 + t197; m(6) * (t192 * t96 - t195 * t97) + m(7) * (t192 * t68 - t195 * t69) + m(5) * t231 * t159; m(7) * (t23 * t24 + t66 * t68 + t67 * t69) + m(6) * (t36 * t47 + t75 * t96 + t76 * t97) + m(5) * (t159 * t209 + t74 * t61) + t216; m(7) * (t24 ^ 2 + t68 ^ 2 + t69 ^ 2) + m(6) * (t47 ^ 2 + t96 ^ 2 + t97 ^ 2) + m(5) * (t159 ^ 2 * t231 + t74 ^ 2) + t216; m(6) * (t62 * t65 + t63 * t64) + m(7) * (t45 * t53 + t46 * t52) + ((t44 / 0.2e1 + t42 / 0.2e1 + t57 / 0.2e1 + t56 / 0.2e1) * t195 + (-t43 / 0.2e1 - t41 / 0.2e1 - t55 / 0.2e1 - t54 / 0.2e1) * t192) * t181 + t261; m(6) * (t192 * t65 - t195 * t64) + m(7) * (t192 * t46 - t195 * t45); m(7) * (t23 * t25 + t45 * t67 + t46 * t66) + m(6) * (t36 * t58 + t64 * t76 + t65 * t75) + t198; m(7) * (t24 * t25 + t45 * t69 + t46 * t68) + m(6) * (t47 * t58 + t64 * t97 + t65 * t96) + t198; t261 * t180 + m(7) * (t25 ^ 2 + t45 ^ 2 + t46 ^ 2) + m(6) * (t58 ^ 2 + t64 ^ 2 + t65 ^ 2) + (t274 * t195 - t275 * t192 + ((t42 + t44) * t195 + (-t41 - t43) * t192) * t180) * t181; m(7) * (t149 * t52 - t151 * t53); m(7) * (t149 * t192 + t151 * t195); m(7) * (t149 * t66 - t151 * t67 + t23 * t247); m(7) * (t149 * t68 - t151 * t69 + t24 * t247); m(7) * (t149 * t46 - t151 * t45 + t247 * t25); m(7) * (t181 ^ 2 * t190 ^ 2 + t149 ^ 2 + t151 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
