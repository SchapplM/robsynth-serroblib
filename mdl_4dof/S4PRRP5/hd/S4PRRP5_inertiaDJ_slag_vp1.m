% Calculate time derivative of joint inertia matrix for
% S4PRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRP5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP5_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_inertiaDJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP5_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP5_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP5_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:28:51
% EndTime: 2019-12-31 16:29:00
% DurationCPUTime: 5.96s
% Computational Cost: add. (4382->381), mult. (12093->584), div. (0->0), fcn. (11481->6), ass. (0->200)
t147 = sin(pkin(6));
t148 = cos(pkin(6));
t152 = cos(qJ(3));
t150 = sin(qJ(3));
t153 = cos(qJ(2));
t213 = t150 * t153;
t138 = t147 * t152 - t148 * t213;
t212 = t152 * t153;
t217 = t147 * t150;
t139 = t148 * t212 + t217;
t151 = sin(qJ(2));
t214 = t148 * t151;
t81 = Icges(5,5) * t139 + Icges(5,6) * t138 + Icges(5,3) * t214;
t83 = Icges(4,5) * t139 + Icges(4,6) * t138 + Icges(4,3) * t214;
t285 = t83 + t81;
t85 = Icges(5,4) * t139 + Icges(5,2) * t138 + Icges(5,6) * t214;
t87 = Icges(4,4) * t139 + Icges(4,2) * t138 + Icges(4,6) * t214;
t284 = t85 + t87;
t89 = Icges(5,1) * t139 + Icges(5,4) * t138 + Icges(5,5) * t214;
t91 = Icges(4,1) * t139 + Icges(4,4) * t138 + Icges(4,5) * t214;
t283 = t89 + t91;
t136 = -t147 * t213 - t148 * t152;
t215 = t148 * t150;
t137 = t147 * t212 - t215;
t216 = t147 * t151;
t84 = Icges(5,4) * t137 + Icges(5,2) * t136 + Icges(5,6) * t216;
t86 = Icges(4,4) * t137 + Icges(4,2) * t136 + Icges(4,6) * t216;
t280 = t84 + t86;
t88 = Icges(5,1) * t137 + Icges(5,4) * t136 + Icges(5,5) * t216;
t90 = Icges(4,1) * t137 + Icges(4,4) * t136 + Icges(4,5) * t216;
t279 = t88 + t90;
t282 = t136 * t284 + t137 * t283 + t216 * t285;
t80 = Icges(5,5) * t137 + Icges(5,6) * t136 + Icges(5,3) * t216;
t82 = Icges(4,5) * t137 + Icges(4,6) * t136 + Icges(4,3) * t216;
t281 = t82 + t80;
t171 = Icges(5,5) * t152 - Icges(5,6) * t150;
t120 = -Icges(5,3) * t153 + t151 * t171;
t172 = Icges(4,5) * t152 - Icges(4,6) * t150;
t121 = -Icges(4,3) * t153 + t151 * t172;
t278 = t120 + t121;
t222 = Icges(5,4) * t152;
t173 = -Icges(5,2) * t150 + t222;
t122 = -Icges(5,6) * t153 + t151 * t173;
t224 = Icges(4,4) * t152;
t174 = -Icges(4,2) * t150 + t224;
t123 = -Icges(4,6) * t153 + t151 * t174;
t274 = t122 + t123;
t223 = Icges(5,4) * t150;
t176 = Icges(5,1) * t152 - t223;
t124 = -Icges(5,5) * t153 + t151 * t176;
t225 = Icges(4,4) * t150;
t177 = Icges(4,1) * t152 - t225;
t125 = -Icges(4,5) * t153 + t151 * t177;
t277 = t124 + t125;
t273 = -t150 * t284 + t152 * t283;
t272 = -t150 * t280 + t152 * t279;
t257 = t136 * t280 + t137 * t279 + t216 * t281;
t271 = t136 * t274 + t137 * t277 + t216 * t278;
t202 = qJD(3) * t151;
t101 = (-Icges(5,2) * t152 - t223) * t202 + (Icges(5,6) * t151 + t153 * t173) * qJD(2);
t102 = (-Icges(4,2) * t152 - t225) * t202 + (Icges(4,6) * t151 + t153 * t174) * qJD(2);
t270 = -t101 - t102;
t103 = (-Icges(5,1) * t150 - t222) * t202 + (Icges(5,5) * t151 + t153 * t176) * qJD(2);
t104 = (-Icges(4,1) * t150 - t224) * t202 + (Icges(4,5) * t151 + t153 * t177) * qJD(2);
t269 = t103 + t104;
t211 = t153 * ((-Icges(4,5) * t150 - Icges(4,6) * t152) * t202 + (Icges(4,3) * t151 + t153 * t172) * qJD(2));
t229 = t153 * ((-Icges(5,5) * t150 - Icges(5,6) * t152) * t202 + (Icges(5,3) * t151 + t153 * t171) * qJD(2));
t266 = -t211 - t229;
t220 = t121 * t153;
t221 = t120 * t153;
t265 = -t220 - t221;
t264 = t282 * t148;
t263 = t151 * t272 - t153 * t281;
t262 = t151 * t273 - t153 * t285;
t146 = t148 ^ 2;
t252 = t147 ^ 2 + t146;
t261 = t150 * t274 - t152 * t277;
t205 = qJD(2) * t151;
t193 = t150 * t205;
t109 = -qJD(3) * t137 + t147 * t193;
t192 = t152 * t205;
t110 = qJD(3) * t136 - t147 * t192;
t204 = qJD(2) * t153;
t195 = t147 * t204;
t57 = Icges(5,5) * t110 + Icges(5,6) * t109 + Icges(5,3) * t195;
t165 = t151 * t57 + t204 * t80;
t61 = Icges(5,4) * t110 + Icges(5,2) * t109 + Icges(5,6) * t195;
t65 = Icges(5,1) * t110 + Icges(5,4) * t109 + Icges(5,5) * t195;
t14 = t109 * t84 + t110 * t88 + t136 * t61 + t137 * t65 + t147 * t165;
t111 = -qJD(3) * t139 + t148 * t193;
t112 = qJD(3) * t138 - t148 * t192;
t194 = t148 * t204;
t58 = Icges(5,5) * t112 + Icges(5,6) * t111 + Icges(5,3) * t194;
t164 = t151 * t58 + t204 * t81;
t62 = Icges(5,4) * t112 + Icges(5,2) * t111 + Icges(5,6) * t194;
t66 = Icges(5,1) * t112 + Icges(5,4) * t111 + Icges(5,5) * t194;
t15 = t109 * t85 + t110 * t89 + t136 * t62 + t137 * t66 + t147 * t164;
t59 = Icges(4,5) * t110 + Icges(4,6) * t109 + Icges(4,3) * t195;
t163 = t151 * t59 + t204 * t82;
t63 = Icges(4,4) * t110 + Icges(4,2) * t109 + Icges(4,6) * t195;
t67 = Icges(4,1) * t110 + Icges(4,4) * t109 + Icges(4,5) * t195;
t16 = t109 * t86 + t110 * t90 + t136 * t63 + t137 * t67 + t147 * t163;
t60 = Icges(4,5) * t112 + Icges(4,6) * t111 + Icges(4,3) * t194;
t162 = t151 * t60 + t204 * t83;
t64 = Icges(4,4) * t112 + Icges(4,2) * t111 + Icges(4,6) * t194;
t68 = Icges(4,1) * t112 + Icges(4,4) * t111 + Icges(4,5) * t194;
t17 = t109 * t87 + t110 * t91 + t136 * t64 + t137 * t68 + t147 * t162;
t260 = (-t109 * t274 - t110 * t277 + t270 * t136 - t269 * t137) * t153 + ((t15 + t17) * t148 + (t14 + t16 + t266) * t147) * t151 + (((t265 + t257) * t147 + t264) * t153 + t271 * t151) * qJD(2);
t259 = (-qJD(2) * t272 + t57 + t59) * t153 + ((-t65 - t67) * t152 + (t61 + t63) * t150 + (t150 * t279 + t152 * t280) * qJD(3) - t281 * qJD(2)) * t151;
t258 = (-qJD(2) * t273 + t58 + t60) * t153 + ((-t66 - t68) * t152 + (t62 + t64) * t150 + (t150 * t283 + t152 * t284) * qJD(3) - t285 * qJD(2)) * t151;
t36 = t138 * t85 + t139 * t89 + t214 * t81;
t38 = t138 * t87 + t139 * t91 + t214 * t83;
t256 = t36 + t38;
t255 = t151 * t261 - t265;
t251 = t147 * t263 + t148 * t262;
t250 = qJD(2) * (rSges(3,1) * t151 + rSges(3,2) * t153);
t239 = pkin(3) * t152;
t248 = qJ(4) * t153 - t151 * t239;
t246 = 2 * m(4);
t245 = 2 * m(5);
t244 = t147 / 0.2e1;
t203 = qJD(3) * t150;
t199 = pkin(3) * t203;
t154 = qJD(2) * t248 + qJD(4) * t151 - t153 * t199;
t201 = qJD(3) * t152;
t198 = pkin(3) * t201;
t237 = rSges(5,1) * t110 + rSges(5,2) * t109 + rSges(5,3) * t195 + t147 * t154 - t148 * t198;
t236 = rSges(5,1) * t112 + rSges(5,2) * t111 + rSges(5,3) * t194 + t147 * t198 + t148 * t154;
t156 = qJ(4) * t151 + t153 * t239;
t235 = rSges(5,1) * t137 + rSges(5,2) * t136 + rSges(5,3) * t216 - pkin(3) * t215 + t147 * t156;
t234 = rSges(5,1) * t139 + rSges(5,2) * t138 + rSges(5,3) * t214 + pkin(3) * t217 + t148 * t156;
t35 = t138 * t84 + t139 * t88 + t214 * t80;
t233 = t147 * t35;
t37 = t138 * t86 + t139 * t90 + t214 * t82;
t232 = t147 * t37;
t186 = rSges(5,1) * t152 - rSges(5,2) * t150;
t228 = -qJD(4) * t153 - t151 * t199 + (-rSges(5,1) * t150 - rSges(5,2) * t152) * t202 + (rSges(5,3) * t151 + t153 * t186 + t156) * qJD(2);
t187 = rSges(4,1) * t152 - rSges(4,2) * t150;
t127 = -rSges(4,3) * t153 + t151 * t187;
t219 = t127 * t153;
t106 = (-rSges(4,1) * t150 - rSges(4,2) * t152) * t202 + (rSges(4,3) * t151 + t153 * t187) * qJD(2);
t189 = pkin(2) * t153 + pkin(5) * t151;
t141 = t189 * qJD(2);
t210 = -t106 - t141;
t209 = -rSges(5,3) * t153 + t151 * t186 - t248;
t143 = pkin(2) * t151 - pkin(5) * t153;
t208 = t252 * qJD(2) * t143;
t207 = -t127 - t143;
t206 = t252 * t189;
t200 = m(5) * t205;
t197 = -t141 - t228;
t196 = -t143 - t209;
t191 = t147 * t209;
t190 = t148 * t209;
t95 = rSges(4,1) * t137 + rSges(4,2) * t136 + rSges(4,3) * t216;
t97 = rSges(4,1) * t139 + rSges(4,2) * t138 + rSges(4,3) * t214;
t183 = -t147 * t97 + t148 * t95;
t157 = qJD(2) * (-Icges(3,5) * t151 - Icges(3,6) * t153);
t155 = -t147 * t234 + t148 * t235;
t130 = t147 * t157;
t108 = t207 * t148;
t107 = t207 * t147;
t79 = t252 * t250;
t78 = t210 * t148;
t77 = t210 * t147;
t76 = t196 * t148;
t75 = t196 * t147;
t72 = rSges(4,1) * t112 + rSges(4,2) * t111 + rSges(4,3) * t194;
t70 = rSges(4,1) * t110 + rSges(4,2) * t109 + rSges(4,3) * t195;
t56 = -t127 * t214 - t153 * t97;
t55 = t127 * t216 + t153 * t95;
t52 = t121 * t214 + t123 * t138 + t125 * t139;
t51 = t120 * t214 + t122 * t138 + t124 * t139;
t48 = t183 * t151;
t47 = t197 * t148;
t46 = t197 * t147;
t45 = t147 * t95 + t148 * t97 + t206;
t40 = -t151 * t190 - t153 * t234;
t39 = t151 * t191 + t153 * t235;
t30 = t147 * t70 + t148 * t72 - t208;
t29 = -t106 * t214 - t153 * t72 + (-t148 * t219 + t151 * t97) * qJD(2);
t28 = t106 * t216 + t153 * t70 + (t147 * t219 - t151 * t95) * qJD(2);
t27 = t155 * t151;
t26 = t147 * t235 + t148 * t234 + t206;
t25 = (-t147 * t72 + t148 * t70) * t151 + t183 * t204;
t24 = t147 * t237 + t148 * t236 - t208;
t23 = -t236 * t153 - t228 * t214 + (t151 * t234 - t153 * t190) * qJD(2);
t22 = t237 * t153 + t228 * t216 + (-t151 * t235 + t153 * t191) * qJD(2);
t21 = t111 * t87 + t112 * t91 + t138 * t64 + t139 * t68 + t148 * t162;
t20 = t111 * t86 + t112 * t90 + t138 * t63 + t139 * t67 + t148 * t163;
t19 = t111 * t85 + t112 * t89 + t138 * t62 + t139 * t66 + t148 * t164;
t18 = t111 * t84 + t112 * t88 + t138 * t61 + t139 * t65 + t148 * t165;
t9 = (-t147 * t236 + t148 * t237) * t151 + t155 * t204;
t8 = t147 * t21 - t148 * t20;
t7 = t147 * t19 - t148 * t18;
t6 = t147 * t17 - t148 * t16;
t5 = -t14 * t148 + t147 * t15;
t4 = -(t138 * t102 + t139 * t104 + t111 * t123 + t112 * t125) * t153 + (t20 * t147 + (t21 - t211) * t148) * t151 + (t52 * t151 + (t232 + (t38 - t220) * t148) * t153) * qJD(2);
t3 = -(t138 * t101 + t139 * t103 + t111 * t122 + t112 * t124) * t153 + (t18 * t147 + (t19 - t229) * t148) * t151 + (t51 * t151 + (t233 + (t36 - t221) * t148) * t153) * qJD(2);
t1 = [0; -m(3) * t79 + m(4) * t30 + m(5) * t24; (t24 * t26 + t46 * t75 + t47 * t76) * t245 + (t107 * t77 + t108 * t78 + t30 * t45) * t246 + 0.2e1 * m(3) * (t250 - t79) * t252 * (rSges(3,1) * t153 - rSges(3,2) * t151) + (-t146 * t130 - t5 - t6) * t148 + (t7 + t8 + (-t147 * t130 + t157 * t252) * t148) * t147; m(4) * t25 + m(5) * t9; t4 * t244 + t3 * t244 + m(5) * (t22 * t76 + t23 * t75 + t24 * t27 + t26 * t9 + t39 * t47 + t40 * t46) + m(4) * (t107 * t29 + t108 * t28 + t25 * t45 + t30 * t48 + t55 * t78 + t56 * t77) + ((t8 / 0.2e1 + t7 / 0.2e1) * t148 + (t6 / 0.2e1 + t5 / 0.2e1) * t147) * t151 + (((t147 * t282 - t257 * t148) * t244 + ((-t35 - t37) * t148 + t256 * t147) * t148 / 0.2e1) * t153 + (t262 * t147 - t263 * t148) * t151 / 0.2e1) * qJD(2) - t260 * t148 / 0.2e1 - (-t258 * t147 + t259 * t148) * t153 / 0.2e1; (t22 * t39 + t23 * t40 + t27 * t9) * t245 + (t25 * t48 + t28 * t55 + t29 * t56) * t246 + t260 * t216 + (t3 + t4) * t214 + (t251 * t205 + (t147 * t257 + t264) * t195 + (t256 * t148 + t232 + t233) * t194) * t151 + (t266 * t153 + t255 * t205 - t271 * t195 + (-t52 - t51) * t194 + ((t150 * t270 + t152 * t269 - t201 * t274 - t203 * t277) * t153 + t258 * t148 + t259 * t147) * t151 + ((-t261 * t153 - t251) * t153 + (t153 * t278 + t255) * t151) * qJD(2)) * t153; t200; m(5) * (-t153 * t24 + (t147 * t46 + t148 * t47) * t151 + (t151 * t26 + (t147 * t75 + t148 * t76) * t153) * qJD(2)); m(5) * (-t153 * t9 + (t147 * t23 + t148 * t22) * t151 + (t151 * t27 + (t147 * t40 + t148 * t39) * t153) * qJD(2)); 0.2e1 * (-0.1e1 + t252) * t153 * t200;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
