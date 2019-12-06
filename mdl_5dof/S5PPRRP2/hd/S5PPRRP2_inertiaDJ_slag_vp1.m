% Calculate time derivative of joint inertia matrix for
% S5PPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:09
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRP2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP2_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP2_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRP2_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRP2_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:08:33
% EndTime: 2019-12-05 15:08:46
% DurationCPUTime: 6.34s
% Computational Cost: add. (8204->389), mult. (12656->594), div. (0->0), fcn. (12148->6), ass. (0->201)
t150 = pkin(8) + qJ(3);
t149 = cos(t150);
t152 = cos(pkin(7));
t153 = sin(qJ(4));
t214 = t152 * t153;
t151 = sin(pkin(7));
t154 = cos(qJ(4));
t215 = t151 * t154;
t138 = t149 * t214 - t215;
t213 = t152 * t154;
t216 = t151 * t153;
t139 = t149 * t213 + t216;
t148 = sin(t150);
t219 = t148 * t152;
t90 = Icges(6,5) * t139 + Icges(6,6) * t219 + Icges(6,3) * t138;
t96 = Icges(5,4) * t139 - Icges(5,2) * t138 + Icges(5,6) * t219;
t285 = t90 - t96;
t92 = Icges(5,5) * t139 - Icges(5,6) * t138 + Icges(5,3) * t219;
t94 = Icges(6,4) * t139 + Icges(6,2) * t219 + Icges(6,6) * t138;
t286 = t92 + t94;
t100 = Icges(5,1) * t139 - Icges(5,4) * t138 + Icges(5,5) * t219;
t98 = Icges(6,1) * t139 + Icges(6,4) * t219 + Icges(6,5) * t138;
t284 = t100 + t98;
t136 = t149 * t216 + t213;
t137 = t149 * t215 - t214;
t220 = t148 * t151;
t89 = Icges(6,5) * t137 + Icges(6,6) * t220 + Icges(6,3) * t136;
t95 = Icges(5,4) * t137 - Icges(5,2) * t136 + Icges(5,6) * t220;
t282 = t89 - t95;
t97 = Icges(6,1) * t137 + Icges(6,4) * t220 + Icges(6,5) * t136;
t99 = Icges(5,1) * t137 - Icges(5,4) * t136 + Icges(5,5) * t220;
t280 = t97 + t99;
t283 = t285 * t136 + t284 * t137 + t286 * t220;
t91 = Icges(5,5) * t137 - Icges(5,6) * t136 + Icges(5,3) * t220;
t93 = Icges(6,4) * t137 + Icges(6,2) * t220 + Icges(6,6) * t136;
t281 = t91 + t93;
t223 = Icges(6,5) * t154;
t172 = Icges(6,3) * t153 + t223;
t113 = -Icges(6,6) * t149 + t148 * t172;
t225 = Icges(5,4) * t154;
t175 = -Icges(5,2) * t153 + t225;
t116 = -Icges(5,6) * t149 + t148 * t175;
t270 = t113 - t116;
t173 = Icges(5,5) * t154 - Icges(5,6) * t153;
t114 = -Icges(5,3) * t149 + t148 * t173;
t174 = Icges(6,4) * t154 + Icges(6,6) * t153;
t115 = -Icges(6,2) * t149 + t148 * t174;
t279 = t114 + t115;
t224 = Icges(6,5) * t153;
t177 = Icges(6,1) * t154 + t224;
t117 = -Icges(6,4) * t149 + t148 * t177;
t226 = Icges(5,4) * t153;
t178 = Icges(5,1) * t154 - t226;
t118 = -Icges(5,5) * t149 + t148 * t178;
t278 = t117 + t118;
t275 = t285 * t153 + t284 * t154;
t274 = t282 * t153 + t280 * t154;
t257 = t282 * t136 + t280 * t137 + t281 * t220;
t273 = t136 * t270 + t278 * t137 + t279 * t220;
t202 = qJD(4) * t148;
t81 = (Icges(6,3) * t154 - t224) * t202 + (Icges(6,6) * t148 + t149 * t172) * qJD(3);
t84 = (-Icges(5,2) * t154 - t226) * t202 + (Icges(5,6) * t148 + t149 * t175) * qJD(3);
t272 = t81 - t84;
t85 = (-Icges(6,1) * t153 + t223) * t202 + (Icges(6,4) * t148 + t149 * t177) * qJD(3);
t86 = (-Icges(5,1) * t153 - t225) * t202 + (Icges(5,5) * t148 + t149 * t178) * qJD(3);
t271 = t85 + t86;
t217 = t149 * t115;
t218 = t149 * t114;
t268 = -t217 - t218;
t267 = t283 * t152;
t234 = t149 * ((-Icges(6,4) * t153 + Icges(6,6) * t154) * t202 + (Icges(6,2) * t148 + t149 * t174) * qJD(3));
t235 = t149 * ((-Icges(5,5) * t153 - Icges(5,6) * t154) * t202 + (Icges(5,3) * t148 + t149 * t173) * qJD(3));
t266 = -t235 - t234;
t265 = t274 * t148 - t149 * t281;
t264 = t275 * t148 - t149 * t286;
t263 = -t270 * t153 - t154 * t278;
t244 = t152 ^ 2;
t251 = t151 ^ 2 + t244;
t200 = qJD(4) * t154;
t193 = t149 * t200;
t204 = qJD(3) * t151;
t109 = -t151 * t193 + (qJD(4) * t152 + t148 * t204) * t153;
t206 = qJD(3) * t148;
t197 = t154 * t206;
t110 = -qJD(4) * t136 - t151 * t197;
t205 = qJD(3) * t149;
t196 = t149 * t204;
t63 = Icges(6,4) * t110 + Icges(6,2) * t196 - Icges(6,6) * t109;
t163 = t148 * t63 + t205 * t93;
t59 = Icges(6,5) * t110 + Icges(6,6) * t196 - Icges(6,3) * t109;
t67 = Icges(6,1) * t110 + Icges(6,4) * t196 - Icges(6,5) * t109;
t14 = -t109 * t89 + t110 * t97 + t136 * t59 + t137 * t67 + t163 * t151;
t201 = qJD(4) * t153;
t203 = qJD(3) * t153;
t111 = -t151 * t201 - t152 * t193 + t203 * t219;
t112 = -qJD(4) * t138 - t152 * t197;
t195 = t152 * t205;
t64 = Icges(6,4) * t112 + Icges(6,2) * t195 - Icges(6,6) * t111;
t162 = t148 * t64 + t205 * t94;
t60 = Icges(6,5) * t112 + Icges(6,6) * t195 - Icges(6,3) * t111;
t68 = Icges(6,1) * t112 + Icges(6,4) * t195 - Icges(6,5) * t111;
t15 = -t109 * t90 + t110 * t98 + t136 * t60 + t137 * t68 + t162 * t151;
t61 = Icges(5,5) * t110 + Icges(5,6) * t109 + Icges(5,3) * t196;
t165 = t148 * t61 + t205 * t91;
t65 = Icges(5,4) * t110 + Icges(5,2) * t109 + Icges(5,6) * t196;
t69 = Icges(5,1) * t110 + Icges(5,4) * t109 + Icges(5,5) * t196;
t16 = t109 * t95 + t110 * t99 - t136 * t65 + t137 * t69 + t165 * t151;
t62 = Icges(5,5) * t112 + Icges(5,6) * t111 + Icges(5,3) * t195;
t164 = t148 * t62 + t205 * t92;
t66 = Icges(5,4) * t112 + Icges(5,2) * t111 + Icges(5,6) * t195;
t70 = Icges(5,1) * t112 + Icges(5,4) * t111 + Icges(5,5) * t195;
t17 = t100 * t110 + t109 * t96 - t136 * t66 + t137 * t70 + t164 * t151;
t262 = (t270 * t109 - t110 * t278 - t272 * t136 - t271 * t137) * t149 + ((t15 + t17) * t152 + (t14 + t16 + t266) * t151) * t148 + (((t268 + t257) * t151 + t267) * t149 + t273 * t148) * qJD(3);
t261 = rSges(6,1) + pkin(4);
t260 = (-t274 * qJD(3) + t61 + t63) * t149 + ((-t67 - t69) * t154 + (-t59 + t65) * t153 + (t280 * t153 - t282 * t154) * qJD(4) - t281 * qJD(3)) * t148;
t259 = (-t275 * qJD(3) + t62 + t64) * t149 + ((-t68 - t70) * t154 + (-t60 + t66) * t153 + (t284 * t153 - t285 * t154) * qJD(4) - t286 * qJD(3)) * t148;
t258 = rSges(6,3) + qJ(5);
t36 = t138 * t90 + t139 * t98 + t94 * t219;
t38 = t100 * t139 - t138 * t96 + t92 * t219;
t256 = t36 + t38;
t255 = t263 * t148 - t268;
t252 = t265 * t151 + t264 * t152;
t250 = qJD(3) * (rSges(4,1) * t148 + rSges(4,2) * t149);
t247 = 2 * m(5);
t246 = 2 * m(6);
t241 = t151 / 0.2e1;
t238 = rSges(6,2) * t196 + qJD(5) * t136 - t258 * t109 + t261 * t110;
t237 = rSges(6,2) * t195 + qJD(5) * t138 - t258 * t111 + t261 * t112;
t186 = pkin(4) * t154 + qJ(5) * t153;
t187 = rSges(6,1) * t154 + rSges(6,3) * t153;
t236 = t186 * t205 + (qJD(5) * t153 + (-pkin(4) * t153 + qJ(5) * t154) * qJD(4)) * t148 + (-rSges(6,1) * t153 + rSges(6,3) * t154) * t202 + (rSges(6,2) * t148 + t149 * t187) * qJD(3);
t35 = t138 * t89 + t139 * t97 + t93 * t219;
t233 = t151 * t35;
t37 = -t138 * t95 + t139 * t99 + t91 * t219;
t232 = t151 * t37;
t190 = pkin(3) * t149 + pkin(6) * t148;
t141 = t190 * qJD(3);
t188 = rSges(5,1) * t154 - rSges(5,2) * t153;
t88 = (-rSges(5,1) * t153 - rSges(5,2) * t154) * t202 + (rSges(5,3) * t148 + t149 * t188) * qJD(3);
t229 = -t141 - t88;
t122 = -t149 * rSges(5,3) + t148 * t188;
t222 = t122 * t149;
t212 = rSges(6,2) * t220 + t258 * t136 + t261 * t137;
t211 = rSges(6,2) * t219 + t258 * t138 + t261 * t139;
t143 = pkin(3) * t148 - pkin(6) * t149;
t210 = t251 * qJD(3) * t143;
t209 = -t149 * rSges(6,2) + (t186 + t187) * t148;
t208 = -t122 - t143;
t207 = t251 * t190;
t199 = -t141 - t236;
t198 = -t143 - t209;
t194 = t149 * t203;
t192 = t209 * t152;
t191 = t209 * t151;
t102 = rSges(5,1) * t137 - rSges(5,2) * t136 + rSges(5,3) * t220;
t104 = rSges(5,1) * t139 - rSges(5,2) * t138 + rSges(5,3) * t219;
t171 = t102 * t152 - t104 * t151;
t157 = qJD(3) * (-Icges(4,5) * t148 - Icges(4,6) * t149);
t156 = t148 * t200 + t194;
t155 = -t211 * t151 + t152 * t212;
t129 = t151 * t157;
t106 = t208 * t152;
t105 = t208 * t151;
t79 = t251 * t250;
t78 = t198 * t152;
t77 = t198 * t151;
t76 = t229 * t152;
t75 = t229 * t151;
t74 = rSges(5,1) * t112 + rSges(5,2) * t111 + rSges(5,3) * t195;
t72 = rSges(5,1) * t110 + rSges(5,2) * t109 + rSges(5,3) * t196;
t58 = -t104 * t149 - t122 * t219;
t57 = t102 * t149 + t122 * t220;
t52 = t171 * t148;
t51 = t114 * t219 - t116 * t138 + t118 * t139;
t50 = t113 * t138 + t115 * t219 + t117 * t139;
t47 = t199 * t152;
t46 = t199 * t151;
t45 = t102 * t151 + t104 * t152 + t207;
t44 = -t148 * t192 - t149 * t211;
t43 = t148 * t191 + t149 * t212;
t30 = t155 * t148;
t29 = t151 * t72 + t152 * t74 - t210;
t28 = -t88 * t219 - t149 * t74 + (t104 * t148 - t152 * t222) * qJD(3);
t27 = t88 * t220 + t149 * t72 + (-t102 * t148 + t151 * t222) * qJD(3);
t26 = t212 * t151 + t152 * t211 + t207;
t25 = (-t151 * t74 + t152 * t72) * t148 + t171 * t205;
t24 = t238 * t151 + t237 * t152 - t210;
t23 = -t237 * t149 - t236 * t219 + (t148 * t211 - t149 * t192) * qJD(3);
t22 = t238 * t149 + t236 * t220 + (-t148 * t212 + t149 * t191) * qJD(3);
t21 = t100 * t112 + t111 * t96 - t138 * t66 + t139 * t70 + t152 * t164;
t20 = t111 * t95 + t112 * t99 - t138 * t65 + t139 * t69 + t152 * t165;
t19 = -t111 * t90 + t112 * t98 + t138 * t60 + t139 * t68 + t152 * t162;
t18 = -t111 * t89 + t112 * t97 + t138 * t59 + t139 * t67 + t152 * t163;
t9 = (-t237 * t151 + t238 * t152) * t148 + t155 * t205;
t8 = t151 * t21 - t152 * t20;
t7 = t151 * t19 - t152 * t18;
t6 = t151 * t17 - t152 * t16;
t5 = -t14 * t152 + t15 * t151;
t4 = -(t111 * t116 + t112 * t118 - t138 * t84 + t139 * t86) * t149 + (t20 * t151 + (t21 - t235) * t152) * t148 + (t51 * t148 + (t232 + (t38 - t218) * t152) * t149) * qJD(3);
t3 = -(-t111 * t113 + t112 * t117 + t138 * t81 + t139 * t85) * t149 + (t18 * t151 + (t19 - t234) * t152) * t148 + (t50 * t148 + (t233 + (t36 - t217) * t152) * t149) * qJD(3);
t1 = [0; 0; 0; -m(4) * t79 + m(5) * t29 + m(6) * t24; m(5) * (t151 * t76 - t152 * t75) + m(6) * (t151 * t47 - t152 * t46); (t26 * t24 + t46 * t77 + t47 * t78) * t246 + (t105 * t75 + t106 * t76 + t45 * t29) * t247 + 0.2e1 * m(4) * (t250 - t79) * t251 * (rSges(4,1) * t149 - rSges(4,2) * t148) + (-t244 * t129 - t5 - t6) * t152 + (t7 + t8 + (-t151 * t129 + t251 * t157) * t152) * t151; m(5) * t25 + m(6) * t9; m(5) * (t151 * t27 - t152 * t28) + m(6) * (t151 * t22 - t152 * t23); t3 * t241 + t4 * t241 + m(6) * (t22 * t78 + t23 * t77 + t30 * t24 + t9 * t26 + t43 * t47 + t44 * t46) + m(5) * (t105 * t28 + t106 * t27 + t25 * t45 + t52 * t29 + t57 * t76 + t58 * t75) + ((t8 / 0.2e1 + t7 / 0.2e1) * t152 + (t6 / 0.2e1 + t5 / 0.2e1) * t151) * t148 + (((t283 * t151 - t257 * t152) * t241 + ((-t35 - t37) * t152 + t256 * t151) * t152 / 0.2e1) * t149 + (t264 * t151 - t265 * t152) * t148 / 0.2e1) * qJD(3) - (-t259 * t151 + t260 * t152) * t149 / 0.2e1 - t262 * t152 / 0.2e1; (t22 * t43 + t23 * t44 + t30 * t9) * t246 + (t25 * t52 + t27 * t57 + t28 * t58) * t247 + t262 * t220 + (t3 + t4) * t219 + (t252 * t206 + (t257 * t151 + t267) * t196 + (t256 * t152 + t232 + t233) * t195) * t148 + (t266 * t149 + t255 * t206 - t273 * t196 + (-t50 - t51) * t195 + ((t272 * t153 + t271 * t154 + t270 * t200 - t201 * t278) * t149 + t259 * t152 + t260 * t151) * t148 + ((-t263 * t149 - t252) * t149 + (t149 * t279 + t255) * t148) * qJD(3)) * t149; t156 * m(6); m(6) * (t109 * t152 - t111 * t151); m(6) * (t26 * t194 - t109 * t77 - t111 * t78 + t136 * t46 + t138 * t47 + (t153 * t24 + t200 * t26) * t148); m(6) * (t30 * t194 - t109 * t44 - t111 * t43 + t136 * t23 + t138 * t22 + (t153 * t9 + t200 * t30) * t148); (t148 * t153 * t156 - t136 * t109 - t138 * t111) * t246;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
