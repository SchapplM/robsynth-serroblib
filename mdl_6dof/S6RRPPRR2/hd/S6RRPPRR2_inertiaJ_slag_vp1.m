% Calculate joint inertia matrix for
% S6RRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
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
% Datum: 2019-03-09 08:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR2_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR2_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRR2_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:49:51
% EndTime: 2019-03-09 08:49:58
% DurationCPUTime: 4.00s
% Computational Cost: add. (9901->464), mult. (8117->676), div. (0->0), fcn. (8597->12), ass. (0->220)
t295 = Icges(3,3) + Icges(4,3);
t196 = qJ(2) + pkin(10);
t186 = sin(t196);
t188 = cos(t196);
t203 = sin(qJ(2));
t205 = cos(qJ(2));
t294 = Icges(3,5) * t205 + Icges(4,5) * t188 - Icges(3,6) * t203 - Icges(4,6) * t186;
t204 = sin(qJ(1));
t282 = t204 / 0.2e1;
t206 = cos(qJ(1));
t293 = t206 / 0.2e1;
t292 = t203 / 0.2e1;
t291 = t205 / 0.2e1;
t200 = cos(pkin(11));
t182 = t200 * pkin(4) + pkin(3);
t195 = pkin(11) + qJ(5);
t187 = cos(t195);
t159 = pkin(5) * t187 + t182;
t185 = sin(t195);
t199 = sin(pkin(11));
t160 = pkin(4) * t199 + pkin(5) * t185;
t261 = t188 * t206;
t189 = qJ(6) + t195;
t180 = sin(t189);
t181 = cos(t189);
t256 = t204 * t181;
t122 = -t180 * t261 + t256;
t257 = t204 * t180;
t123 = t181 * t261 + t257;
t262 = t186 * t206;
t73 = t123 * rSges(7,1) + t122 * rSges(7,2) + rSges(7,3) * t262;
t290 = t159 * t261 + t204 * t160 + t73;
t289 = -t294 * t204 + t295 * t206;
t288 = t295 * t204 + t294 * t206;
t197 = t204 ^ 2;
t198 = t206 ^ 2;
t287 = m(5) / 0.2e1;
t286 = m(6) / 0.2e1;
t285 = m(7) / 0.2e1;
t263 = t186 * t204;
t120 = -t181 * t206 - t188 * t257;
t121 = -t180 * t206 + t188 * t256;
t66 = Icges(7,5) * t121 + Icges(7,6) * t120 + Icges(7,3) * t263;
t68 = Icges(7,4) * t121 + Icges(7,2) * t120 + Icges(7,6) * t263;
t70 = Icges(7,1) * t121 + Icges(7,4) * t120 + Icges(7,5) * t263;
t22 = t120 * t68 + t121 * t70 + t263 * t66;
t67 = Icges(7,5) * t123 + Icges(7,6) * t122 + Icges(7,3) * t262;
t69 = Icges(7,4) * t123 + Icges(7,2) * t122 + Icges(7,6) * t262;
t71 = Icges(7,1) * t123 + Icges(7,4) * t122 + Icges(7,5) * t262;
t23 = t120 * t69 + t121 * t71 + t263 * t67;
t100 = -Icges(7,3) * t188 + (Icges(7,5) * t181 - Icges(7,6) * t180) * t186;
t101 = -Icges(7,6) * t188 + (Icges(7,4) * t181 - Icges(7,2) * t180) * t186;
t102 = -Icges(7,5) * t188 + (Icges(7,1) * t181 - Icges(7,4) * t180) * t186;
t39 = t100 * t263 + t101 * t120 + t102 * t121;
t5 = -t39 * t188 + (t204 * t22 + t206 * t23) * t186;
t24 = t122 * t68 + t123 * t70 + t262 * t66;
t25 = t122 * t69 + t123 * t71 + t262 * t67;
t40 = t100 * t262 + t122 * t101 + t123 * t102;
t6 = -t40 * t188 + (t204 * t24 + t206 * t25) * t186;
t284 = t6 * t262 + t5 * t263;
t283 = -t188 / 0.2e1;
t281 = -t206 / 0.2e1;
t166 = rSges(3,1) * t203 + rSges(3,2) * t205;
t280 = m(3) * t166;
t279 = pkin(2) * t203;
t278 = pkin(3) * t188;
t277 = -pkin(3) + t182;
t202 = -pkin(8) - qJ(4);
t194 = -pkin(9) + t202;
t244 = t194 - t202;
t253 = t204 * t199;
t249 = -pkin(4) * t253 - t182 * t261;
t276 = -t244 * t262 + t249 + t290;
t103 = -rSges(7,3) * t188 + (rSges(7,1) * t181 - rSges(7,2) * t180) * t186;
t222 = -t121 * rSges(7,1) - t120 * rSges(7,2);
t72 = rSges(7,3) * t263 - t222;
t52 = t103 * t263 + t188 * t72;
t275 = rSges(3,1) * t205;
t274 = rSges(3,2) * t203;
t273 = t206 * rSges(3,3);
t265 = t101 * t180;
t94 = t186 * t181 * t102;
t47 = -t188 * t100 - t186 * t265 + t94;
t272 = t47 * t188;
t248 = t159 - t182;
t92 = t186 * t248 + t188 * t244;
t271 = -t103 - t92;
t270 = Icges(3,4) * t203;
t269 = Icges(3,4) * t205;
t268 = Icges(4,4) * t186;
t267 = Icges(4,4) * t188;
t266 = qJ(4) * t186;
t105 = -Icges(6,6) * t188 + (Icges(6,4) * t187 - Icges(6,2) * t185) * t186;
t264 = t105 * t185;
t260 = t199 * t206;
t259 = t200 * t206;
t201 = -qJ(3) - pkin(7);
t258 = t201 * t206;
t255 = t204 * t185;
t254 = t204 * t187;
t252 = t204 * t200;
t183 = pkin(2) * t205 + pkin(1);
t175 = t206 * t183;
t193 = t206 * pkin(7);
t251 = t204 * (t258 + t193 + (-pkin(1) + t183) * t204) + t206 * (-pkin(1) * t206 + t175 + (-pkin(7) - t201) * t204);
t247 = pkin(4) * t260 + t202 * t263;
t246 = pkin(3) * t261 + qJ(4) * t262;
t245 = t204 * rSges(3,3) + t206 * t275;
t243 = t197 + t198;
t135 = -t187 * t206 - t188 * t255;
t136 = -t185 * t206 + t188 * t254;
t74 = Icges(6,5) * t136 + Icges(6,6) * t135 + Icges(6,3) * t263;
t76 = Icges(6,4) * t136 + Icges(6,2) * t135 + Icges(6,6) * t263;
t78 = Icges(6,1) * t136 + Icges(6,4) * t135 + Icges(6,5) * t263;
t35 = -t188 * t74 + (-t185 * t76 + t187 * t78) * t186;
t104 = -Icges(6,3) * t188 + (Icges(6,5) * t187 - Icges(6,6) * t185) * t186;
t106 = -Icges(6,5) * t188 + (Icges(6,1) * t187 - Icges(6,4) * t185) * t186;
t43 = t104 * t263 + t105 * t135 + t106 * t136;
t242 = t35 / 0.2e1 + t43 / 0.2e1;
t137 = -t185 * t261 + t254;
t138 = t187 * t261 + t255;
t75 = Icges(6,5) * t138 + Icges(6,6) * t137 + Icges(6,3) * t262;
t77 = Icges(6,4) * t138 + Icges(6,2) * t137 + Icges(6,6) * t262;
t79 = Icges(6,1) * t138 + Icges(6,4) * t137 + Icges(6,5) * t262;
t36 = -t188 * t75 + (-t185 * t77 + t187 * t79) * t186;
t44 = t104 * t262 + t137 * t105 + t138 * t106;
t241 = t36 / 0.2e1 + t44 / 0.2e1;
t81 = t138 * rSges(6,1) + t137 * rSges(6,2) + rSges(6,3) * t262;
t150 = -t188 * t260 + t252;
t151 = t188 * t259 + t253;
t240 = t151 * rSges(5,1) + t150 * rSges(5,2) + rSges(5,3) * t262;
t239 = t263 / 0.2e1;
t238 = t262 / 0.2e1;
t237 = Icges(4,5) * t186 / 0.2e1 + Icges(4,6) * t188 / 0.2e1 + Icges(3,5) * t292 + Icges(3,6) * t291;
t236 = -pkin(3) * t186 + qJ(4) * t188 - t279;
t235 = -rSges(4,1) * t186 - rSges(4,2) * t188 - t279;
t32 = -t188 * t66 + (-t180 * t68 + t181 * t70) * t186;
t33 = -t188 * t67 + (-t180 * t69 + t181 * t71) * t186;
t234 = (t32 + t39) * t239 + (t33 + t40) * t238;
t233 = -t204 * t201 + t175;
t9 = -t272 + (t204 * t32 + t206 * t33) * t186;
t232 = -t188 * t9 + t284;
t231 = t287 + t286 + t285;
t230 = t197 * (t266 + t278) + t206 * t246 + t251;
t13 = t23 * t204 - t206 * t22;
t14 = t25 * t204 - t206 * t24;
t229 = t13 * t239 + t14 * t238 + t5 * t281 + t6 * t282 + (t33 * t204 - t32 * t206) * t283;
t228 = t236 - (qJ(4) + t202) * t188 - t277 * t186;
t227 = rSges(5,3) * t188 - (rSges(5,1) * t200 - rSges(5,2) * t199) * t186 + t236;
t226 = -t274 + t275;
t225 = rSges(4,1) * t188 - rSges(4,2) * t186;
t148 = -t188 * t253 - t259;
t149 = t188 * t252 - t260;
t224 = -t149 * rSges(5,1) - t148 * rSges(5,2);
t223 = -t136 * rSges(6,1) - t135 * rSges(6,2);
t221 = Icges(3,1) * t205 - t270;
t220 = Icges(4,1) * t188 - t268;
t219 = -Icges(3,2) * t203 + t269;
t218 = -Icges(4,2) * t186 + t267;
t107 = -rSges(6,3) * t188 + (rSges(6,1) * t187 - rSges(6,2) * t185) * t186;
t211 = -t107 + t228;
t210 = rSges(4,1) * t261 - rSges(4,2) * t262 + t204 * rSges(4,3);
t209 = -t202 * t262 - t249;
t208 = t204 * ((t188 * t277 - t266) * t204 - t247) + t206 * (t209 - t246) + t230;
t207 = t228 + t271;
t168 = rSges(2,1) * t206 - t204 * rSges(2,2);
t167 = -t204 * rSges(2,1) - rSges(2,2) * t206;
t119 = t235 * t206;
t118 = t235 * t204;
t114 = -Icges(5,5) * t188 + (Icges(5,1) * t200 - Icges(5,4) * t199) * t186;
t113 = -Icges(5,6) * t188 + (Icges(5,4) * t200 - Icges(5,2) * t199) * t186;
t111 = t204 * pkin(7) + (pkin(1) - t274) * t206 + t245;
t110 = t273 + t193 + (-pkin(1) - t226) * t204;
t98 = t210 + t233;
t97 = (rSges(4,3) - t201) * t206 + (-t183 - t225) * t204;
t95 = t186 * t187 * t106;
t93 = t206 * (-t206 * t274 + t245) + (t204 * t226 - t273) * t204;
t91 = Icges(5,1) * t151 + Icges(5,4) * t150 + Icges(5,5) * t262;
t90 = Icges(5,1) * t149 + Icges(5,4) * t148 + Icges(5,5) * t263;
t89 = Icges(5,4) * t151 + Icges(5,2) * t150 + Icges(5,6) * t262;
t88 = Icges(5,4) * t149 + Icges(5,2) * t148 + Icges(5,6) * t263;
t87 = Icges(5,5) * t151 + Icges(5,6) * t150 + Icges(5,3) * t262;
t86 = Icges(5,5) * t149 + Icges(5,6) * t148 + Icges(5,3) * t263;
t85 = t227 * t206;
t84 = t227 * t204;
t80 = rSges(6,3) * t263 - t223;
t64 = t72 * t262;
t62 = -t160 * t206 + (-t186 * t194 + t188 * t248) * t204 + t247;
t61 = t233 + t240 + t246;
t60 = -t258 + (-t278 - t183 + (-rSges(5,3) - qJ(4)) * t186) * t204 + t224;
t59 = t211 * t206;
t58 = t211 * t204;
t57 = -t107 * t262 - t188 * t81;
t56 = t107 * t263 + t188 * t80;
t55 = t209 + t233 + t81;
t54 = -t258 + (-rSges(6,3) * t186 - t182 * t188 - t183) * t204 + t223 + t247;
t53 = -t103 * t262 - t188 * t73;
t51 = t206 * t210 + (-t206 * rSges(4,3) + t204 * t225) * t204 + t251;
t50 = -t188 * t104 - t186 * t264 + t95;
t49 = -t194 * t262 + t233 + t290;
t48 = (t160 - t201) * t206 + (-t159 * t188 - t183 + (-rSges(7,3) + t194) * t186) * t204 + t222;
t46 = (-t204 * t81 + t206 * t80) * t186;
t45 = -t263 * t73 + t64;
t42 = t207 * t206;
t41 = t207 * t204;
t34 = t204 * (rSges(5,3) * t263 - t224) + t206 * t240 + t230;
t29 = t137 * t77 + t138 * t79 + t262 * t75;
t28 = t137 * t76 + t138 * t78 + t262 * t74;
t27 = t135 * t77 + t136 * t79 + t263 * t75;
t26 = t135 * t76 + t136 * t78 + t263 * t74;
t21 = -t188 * t276 + t262 * t271;
t20 = t188 * t62 + t263 * t92 + t52;
t19 = t64 + (-t204 * t276 + t206 * t62) * t186;
t18 = t204 * t80 + t206 * t81 + t208;
t16 = t29 * t204 - t206 * t28;
t15 = t27 * t204 - t206 * t26;
t12 = t276 * t206 + (t62 + t72) * t204 + t208;
t8 = -t44 * t188 + (t204 * t28 + t206 * t29) * t186;
t7 = -t43 * t188 + (t204 * t26 + t206 * t27) * t186;
t1 = [t205 * (Icges(3,2) * t205 + t270) + t203 * (Icges(3,1) * t203 + t269) + Icges(2,3) + t94 + t95 + (-t100 - t104 + t268 - (Icges(5,5) * t200 - Icges(5,6) * t199) * t186 + (Icges(4,2) + Icges(5,3)) * t188) * t188 + (Icges(4,1) * t186 - t113 * t199 + t114 * t200 - t264 - t265 + t267) * t186 + m(7) * (t48 ^ 2 + t49 ^ 2) + m(6) * (t54 ^ 2 + t55 ^ 2) + m(5) * (t60 ^ 2 + t61 ^ 2) + m(4) * (t97 ^ 2 + t98 ^ 2) + m(3) * (t110 ^ 2 + t111 ^ 2) + m(2) * (t167 ^ 2 + t168 ^ 2); m(7) * (t41 * t49 + t42 * t48) + m(6) * (t54 * t59 + t55 * t58) + m(5) * (t60 * t85 + t61 * t84) + m(4) * (t118 * t98 + t119 * t97) + (-t113 * t148 / 0.2e1 - t114 * t149 / 0.2e1 - t110 * t280 - t32 / 0.2e1 - t39 / 0.2e1 - (-Icges(3,6) * t206 + t204 * t219) * t205 / 0.2e1 - (-Icges(3,5) * t206 + t204 * t221) * t203 / 0.2e1 + t237 * t206 + (Icges(4,6) * t293 - t204 * t218 / 0.2e1 + t86 / 0.2e1) * t188 - t242) * t206 + (t150 * t113 / 0.2e1 + t151 * t114 / 0.2e1 - t111 * t280 + t33 / 0.2e1 + t40 / 0.2e1 + (Icges(3,6) * t204 + t206 * t219) * t291 + (Icges(3,5) * t204 + t206 * t221) * t292 + t237 * t204 + (Icges(4,6) * t282 + t218 * t293 - t87 / 0.2e1) * t188 + t241) * t204 + ((Icges(4,5) * t204 - t199 * t89 + t200 * t91 + t206 * t220) * t282 + (-Icges(4,5) * t206 - t199 * t88 + t200 * t90 + t204 * t220) * t281) * t186; m(7) * (t12 ^ 2 + t41 ^ 2 + t42 ^ 2) + m(6) * (t18 ^ 2 + t58 ^ 2 + t59 ^ 2) + m(5) * (t34 ^ 2 + t84 ^ 2 + t85 ^ 2) + m(4) * (t118 ^ 2 + t119 ^ 2 + t51 ^ 2) + m(3) * (t166 ^ 2 * t243 + t93 ^ 2) + (-t13 - t15 + (t148 * t88 + t149 * t90 + t86 * t263) * t206 + t289 * t198) * t206 + (t14 + t16 + (t150 * t89 + t151 * t91 + t87 * t262) * t204 + t288 * t197 + (-t148 * t89 - t149 * t91 - t150 * t88 - t151 * t90 + t204 * t289 + t206 * t288 - t86 * t262 - t87 * t263) * t206) * t204; m(7) * (t204 * t48 - t206 * t49) + m(6) * (t204 * t54 - t206 * t55) + m(5) * (t204 * t60 - t206 * t61) + m(4) * (t204 * t97 - t206 * t98); m(7) * (t204 * t42 - t206 * t41) + m(6) * (t204 * t59 - t206 * t58) + m(5) * (t204 * t85 - t206 * t84) + m(4) * (-t118 * t206 + t204 * t119); 0.2e1 * (m(4) / 0.2e1 + t231) * t243; 0.2e1 * ((t204 * t49 + t206 * t48) * t285 + (t204 * t55 + t206 * t54) * t286 + (t204 * t61 + t206 * t60) * t287) * t186; m(7) * (-t188 * t12 + (t204 * t41 + t206 * t42) * t186) + m(6) * (-t188 * t18 + (t204 * t58 + t206 * t59) * t186) + m(5) * (-t188 * t34 + (t204 * t84 + t206 * t85) * t186); 0; 0.2e1 * t231 * (t186 ^ 2 * t243 + t188 ^ 2); (-t47 - t50) * t188 + m(7) * (t20 * t48 + t21 * t49) + m(6) * (t54 * t56 + t55 * t57) + (t204 * t242 + t206 * t241) * t186 + t234; (t36 * t204 - t35 * t206) * t283 + t7 * t281 + t8 * t282 + (t15 * t282 + t16 * t293) * t186 + m(7) * (t12 * t19 + t20 * t42 + t21 * t41) + m(6) * (t18 * t46 + t56 * t59 + t57 * t58) + t229; m(6) * (t56 * t204 - t206 * t57) + m(7) * (t20 * t204 - t206 * t21); m(6) * (-t46 * t188 + (t204 * t57 + t206 * t56) * t186) + m(7) * (-t19 * t188 + (t20 * t206 + t204 * t21) * t186); (t50 * t188 - t9) * t188 + m(7) * (t19 ^ 2 + t20 ^ 2 + t21 ^ 2) + m(6) * (t46 ^ 2 + t56 ^ 2 + t57 ^ 2) + (t206 * t8 + t204 * t7 - t188 * (t204 * t35 + t206 * t36)) * t186 + t284; -t272 + m(7) * (t48 * t52 + t49 * t53) + t234; m(7) * (t12 * t45 + t41 * t53 + t42 * t52) + t229; m(7) * (t52 * t204 - t206 * t53); m(7) * (-t45 * t188 + (t204 * t53 + t206 * t52) * t186); m(7) * (t19 * t45 + t20 * t52 + t21 * t53) + t232; m(7) * (t45 ^ 2 + t52 ^ 2 + t53 ^ 2) + t232;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
