% Calculate time derivative of joint inertia matrix for
% S5RPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-12-31 17:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR6_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR6_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR6_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR6_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR6_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR6_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR6_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:57:41
% EndTime: 2019-12-31 17:57:49
% DurationCPUTime: 4.12s
% Computational Cost: add. (13102->413), mult. (12677->623), div. (0->0), fcn. (12061->10), ass. (0->216)
t141 = pkin(9) + qJ(4);
t136 = sin(t141);
t138 = cos(t141);
t146 = sin(qJ(5));
t148 = cos(qJ(5));
t228 = Icges(6,4) * t148;
t167 = -Icges(6,2) * t146 + t228;
t101 = -Icges(6,6) * t138 + t167 * t136;
t229 = Icges(6,4) * t146;
t170 = Icges(6,1) * t148 - t229;
t102 = -Icges(6,5) * t138 + t170 * t136;
t263 = -t101 * t148 - t102 * t146;
t165 = Icges(6,5) * t148 - Icges(6,6) * t146;
t100 = -Icges(6,3) * t138 + t165 * t136;
t205 = qJD(4) * t148;
t206 = qJD(4) * t146;
t210 = qJD(4) * t136;
t204 = qJD(5) * t136;
t82 = (-Icges(6,5) * t146 - Icges(6,6) * t148) * t204 + (Icges(6,3) * t136 + t165 * t138) * qJD(4);
t84 = (-Icges(6,1) * t146 - t228) * t204 + (Icges(6,5) * t136 + t170 * t138) * qJD(4);
t262 = t136 * t148 * t84 + t100 * t210 + (-t101 * t206 + t102 * t205 - t82) * t138;
t142 = qJ(1) + pkin(8);
t137 = sin(t142);
t139 = cos(t142);
t230 = Icges(5,4) * t138;
t169 = -Icges(5,2) * t136 + t230;
t95 = Icges(5,6) * t137 + t169 * t139;
t231 = Icges(5,4) * t136;
t172 = Icges(5,1) * t138 - t231;
t97 = Icges(5,5) * t137 + t172 * t139;
t180 = t136 * t95 - t138 * t97;
t261 = t180 * t137;
t94 = -Icges(5,6) * t139 + t169 * t137;
t96 = -Icges(5,5) * t139 + t172 * t137;
t181 = t136 * t94 - t138 * t96;
t260 = t181 * t139;
t134 = t137 * rSges(5,3);
t220 = t136 * t139;
t259 = -rSges(5,2) * t220 + t134;
t145 = -pkin(6) - qJ(3);
t144 = cos(pkin(9));
t135 = pkin(3) * t144 + pkin(2);
t243 = rSges(5,1) * t138;
t187 = -rSges(5,2) * t136 + t243;
t162 = -t135 - t187;
t247 = sin(qJ(1)) * pkin(1);
t80 = -t247 + (rSges(5,3) - t145) * t139 + t162 * t137;
t140 = cos(qJ(1)) * pkin(1);
t189 = t139 * t135 - t137 * t145 + t140;
t217 = t138 * t139;
t99 = rSges(5,1) * t217 + t259;
t81 = t189 + t99;
t258 = t137 * t81 + t139 * t80;
t166 = Icges(5,5) * t138 - Icges(5,6) * t136;
t92 = -Icges(5,3) * t139 + t166 * t137;
t163 = rSges(4,1) * t144 - rSges(4,2) * sin(pkin(9)) + pkin(2);
t235 = rSges(4,3) + qJ(3);
t89 = t235 * t137 + t163 * t139 + t140;
t257 = 2 * m(5);
t256 = 2 * m(6);
t255 = t137 ^ 2;
t254 = t139 ^ 2;
t253 = t137 / 0.2e1;
t252 = -t138 / 0.2e1;
t250 = t139 / 0.2e1;
t249 = -rSges(6,3) - pkin(7);
t120 = rSges(5,1) * t136 + rSges(5,2) * t138;
t248 = m(5) * t120;
t246 = pkin(4) * t136;
t245 = pkin(4) * t138;
t244 = qJD(1) / 0.2e1;
t242 = t136 * t96;
t241 = t136 * t97;
t239 = t138 * t94;
t238 = t138 * t95;
t83 = (-Icges(6,2) * t148 - t229) * t204 + (Icges(6,6) * t136 + t167 * t138) * qJD(4);
t236 = t146 * t83;
t190 = pkin(7) * t136 + t245;
t215 = t139 * t148;
t219 = t137 * t146;
t108 = -t138 * t219 - t215;
t216 = t139 * t146;
t218 = t137 * t148;
t109 = t138 * t218 - t216;
t186 = -rSges(6,1) * t109 - rSges(6,2) * t108;
t221 = t136 * t137;
t76 = rSges(6,3) * t221 - t186;
t234 = t190 * t137 + t76;
t110 = -t138 * t216 + t218;
t111 = t138 * t215 + t219;
t77 = t111 * rSges(6,1) + t110 * rSges(6,2) + rSges(6,3) * t220;
t233 = pkin(4) * t217 + pkin(7) * t220 + t77;
t185 = rSges(6,1) * t148 - rSges(6,2) * t146;
t85 = (-rSges(6,1) * t146 - rSges(6,2) * t148) * t204 + (rSges(6,3) * t136 + t185 * t138) * qJD(4);
t232 = -t190 * qJD(4) - t85;
t93 = Icges(5,3) * t137 + t166 * t139;
t224 = qJD(1) * t93;
t103 = -rSges(6,3) * t138 + t185 * t136;
t121 = -pkin(7) * t138 + t246;
t214 = t103 + t121;
t133 = qJD(3) * t139;
t212 = qJD(1) * t137;
t213 = t145 * t212 + t133;
t211 = qJD(1) * t139;
t209 = qJD(4) * t137;
t208 = qJD(4) * t138;
t207 = qJD(4) * t139;
t195 = t138 * t207;
t193 = -qJD(5) * t138 + qJD(1);
t153 = t136 * t206 + t193 * t148;
t192 = t138 * qJD(1) - qJD(5);
t59 = t153 * t139 + t192 * t219;
t152 = -t136 * t205 + t193 * t146;
t60 = t152 * t139 - t192 * t218;
t202 = t60 * rSges(6,1) + t59 * rSges(6,2) + rSges(6,3) * t195;
t72 = Icges(6,4) * t109 + Icges(6,2) * t108 + Icges(6,6) * t221;
t74 = Icges(6,1) * t109 + Icges(6,4) * t108 + Icges(6,5) * t221;
t174 = -t146 * t72 + t148 * t74;
t70 = Icges(6,5) * t109 + Icges(6,6) * t108 + Icges(6,3) * t221;
t29 = t174 * t136 - t138 * t70;
t41 = t100 * t221 + t101 * t108 + t102 * t109;
t200 = t29 / 0.2e1 + t41 / 0.2e1;
t73 = Icges(6,4) * t111 + Icges(6,2) * t110 + Icges(6,6) * t220;
t75 = Icges(6,1) * t111 + Icges(6,4) * t110 + Icges(6,5) * t220;
t173 = -t146 * t73 + t148 * t75;
t71 = Icges(6,5) * t111 + Icges(6,6) * t110 + Icges(6,3) * t220;
t30 = t173 * t136 - t138 * t71;
t42 = t100 * t220 + t101 * t110 + t102 * t111;
t199 = t30 / 0.2e1 + t42 / 0.2e1;
t198 = t136 * t212;
t196 = t137 * t208;
t194 = t208 / 0.2e1;
t87 = t214 * t139;
t61 = t153 * t137 - t192 * t216;
t62 = t152 * t137 + t192 * t215;
t191 = t62 * rSges(6,1) + t61 * rSges(6,2);
t184 = -t139 * t145 - t247;
t24 = t108 * t72 + t109 * t74 + t70 * t221;
t25 = t108 * t73 + t109 * t75 + t71 * t221;
t15 = t137 * t25 - t139 * t24;
t179 = t137 * t24 + t139 * t25;
t26 = t110 * t72 + t111 * t74 + t70 * t220;
t27 = t110 * t73 + t111 * t75 + t71 * t220;
t16 = t137 * t27 - t139 * t26;
t178 = t137 * t26 + t139 * t27;
t177 = t30 * t137 - t29 * t139;
t176 = t137 * t29 + t139 * t30;
t175 = -t137 * t77 + t139 * t76;
t171 = Icges(5,1) * t136 + t230;
t168 = Icges(5,2) * t138 + t231;
t164 = pkin(7) * t195 - t207 * t246;
t161 = qJD(4) * t120;
t160 = qJD(4) * t171;
t159 = qJD(4) * t168;
t158 = qJD(4) * (-Icges(5,5) * t136 - Icges(5,6) * t138);
t157 = t249 * t136 - t135 - t245;
t155 = t136 * t211 + t196;
t154 = t195 - t198;
t151 = rSges(5,2) * t198 + rSges(5,3) * t211 - t139 * t161;
t150 = t157 * t137 + t184;
t88 = -t163 * t137 + t235 * t139 - t247;
t132 = qJD(3) * t137;
t115 = t187 * qJD(4);
t98 = -rSges(5,3) * t139 + t187 * t137;
t86 = t214 * t137;
t79 = -t89 * qJD(1) + t133;
t78 = t88 * qJD(1) + t132;
t64 = t137 * t158 + t224;
t63 = -qJD(1) * t92 + t139 * t158;
t56 = t120 * t209 + (t162 * t139 - t134 - t140) * qJD(1) + t213;
t55 = t132 + ((-t135 - t243) * t137 + t184) * qJD(1) + t151;
t54 = -t100 * t138 + (-t101 * t146 + t102 * t148) * t136;
t53 = -t103 * t220 - t138 * t77;
t52 = t103 * t221 + t138 * t76;
t51 = t189 + t233;
t50 = t150 + t186;
t49 = t54 * t210;
t48 = t137 * t93 - t180 * t139;
t47 = t137 * t92 - t260;
t46 = -t139 * t93 - t261;
t45 = -t181 * t137 - t139 * t92;
t44 = -qJD(1) * t87 + t232 * t137;
t43 = t232 * t139 + t214 * t212;
t40 = t175 * t136;
t39 = t155 * rSges(6,3) + t191;
t38 = -rSges(6,3) * t198 + t202;
t37 = Icges(6,1) * t62 + Icges(6,4) * t61 + t155 * Icges(6,5);
t36 = Icges(6,1) * t60 + Icges(6,4) * t59 + t154 * Icges(6,5);
t35 = Icges(6,4) * t62 + Icges(6,2) * t61 + t155 * Icges(6,6);
t34 = Icges(6,4) * t60 + Icges(6,2) * t59 + t154 * Icges(6,6);
t33 = Icges(6,5) * t62 + Icges(6,6) * t61 + t155 * Icges(6,3);
t32 = Icges(6,5) * t60 + Icges(6,6) * t59 + t154 * Icges(6,3);
t31 = t234 * t137 + t233 * t139;
t28 = (qJD(1) * t98 + t151) * t139 + (-t137 * t161 + (-t99 + t259) * qJD(1)) * t137;
t23 = (t249 * t138 + t246) * t209 + (t157 * t139 - t140) * qJD(1) - t191 + t213;
t22 = t150 * qJD(1) + t132 + t164 + t202;
t21 = (t263 * qJD(5) - t236) * t136 + t262;
t20 = (t103 * t209 + t39) * t138 + (-qJD(4) * t76 + t103 * t211 + t137 * t85) * t136;
t19 = (-t103 * t207 - t38) * t138 + (qJD(4) * t77 + t103 * t212 - t139 * t85) * t136;
t18 = t155 * t100 + t101 * t61 + t102 * t62 + t108 * t83 + t109 * t84 + t82 * t221;
t17 = t154 * t100 + t101 * t59 + t102 * t60 + t110 * t83 + t111 * t84 + t82 * t220;
t14 = t175 * t208 + (-t137 * t38 + t139 * t39 + (-t137 * t76 - t139 * t77) * qJD(1)) * t136;
t13 = (t164 + t38) * t139 + (-t121 * t209 + t39) * t137 + (-t233 * t137 + t234 * t139) * qJD(1);
t12 = t178 * t136 - t138 * t42;
t11 = t179 * t136 - t138 * t41;
t10 = (t173 * qJD(4) - t32) * t138 + (qJD(4) * t71 - t146 * t34 + t148 * t36 + (-t146 * t75 - t148 * t73) * qJD(5)) * t136;
t9 = (t174 * qJD(4) - t33) * t138 + (qJD(4) * t70 - t146 * t35 + t148 * t37 + (-t146 * t74 - t148 * t72) * qJD(5)) * t136;
t8 = t71 * t196 + t108 * t34 + t109 * t36 + t61 * t73 + t62 * t75 + (t137 * t32 + t71 * t211) * t136;
t7 = t70 * t196 + t108 * t35 + t109 * t37 + t61 * t72 + t62 * t74 + (t137 * t33 + t70 * t211) * t136;
t6 = t71 * t195 + t110 * t34 + t111 * t36 + t59 * t73 + t60 * t75 + (t139 * t32 - t71 * t212) * t136;
t5 = t70 * t195 + t110 * t35 + t111 * t37 + t59 * t72 + t60 * t74 + (t139 * t33 - t70 * t212) * t136;
t4 = t179 * qJD(1) + t137 * t8 - t139 * t7;
t3 = t178 * qJD(1) + t137 * t6 - t139 * t5;
t2 = (t179 * qJD(4) - t18) * t138 + (-t15 * qJD(1) + qJD(4) * t41 + t137 * t7 + t139 * t8) * t136;
t1 = (t178 * qJD(4) - t17) * t138 + (-t16 * qJD(1) + qJD(4) * t42 + t137 * t5 + t139 * t6) * t136;
t57 = [(t22 * t51 + t23 * t50) * t256 - t136 * t236 + (t55 * t81 + t56 * t80) * t257 + 0.2e1 * m(4) * (t78 * t89 + t79 * t88) + (t172 - t168) * t210 + (t171 + t169) * t208 + t263 * t204 + t262; 0; 0; m(6) * (t137 * t23 - t139 * t22 + (t137 * t51 + t139 * t50) * qJD(1)) + m(5) * (t258 * qJD(1) + t137 * t56 - t139 * t55) + m(4) * (t137 * t79 - t139 * t78 + (t137 * t89 + t139 * t88) * qJD(1)); 0; 0; m(6) * (-t22 * t86 - t23 * t87 + t43 * t50 + t44 * t51) + m(5) * ((-t137 * t55 - t139 * t56) * t120 - t258 * t115) + (t255 / 0.2e1 + t254 / 0.2e1) * t166 * qJD(4) + ((t241 / 0.2e1 + t238 / 0.2e1 - t81 * t248 + t199) * t139 + (t242 / 0.2e1 + t239 / 0.2e1 + t80 * t248 + t200) * t137) * qJD(1) + (-t180 * qJD(4) + t136 * (-t96 * qJD(1) - t139 * t160) + t138 * (-t94 * qJD(1) - t139 * t159) + t10 + t17) * t253 - (-t181 * qJD(4) + t136 * (t97 * qJD(1) - t137 * t160) + t138 * (t95 * qJD(1) - t137 * t159) + t18 + t9) * t139 / 0.2e1; m(5) * t28 + m(6) * t13; m(6) * (t137 * t43 - t139 * t44 + (-t137 * t86 - t139 * t87) * qJD(1)); ((t137 * t98 + t139 * t99) * t28 + (t254 + t255) * t120 * t115) * t257 + t137 * ((t137 * t63 + (t47 + t261) * qJD(1)) * t137 + (t48 * qJD(1) + (t94 * t208 + t96 * t210) * t139 + (-t64 + (-t238 - t241) * qJD(4) + (-t181 + t93) * qJD(1)) * t137) * t139) - t139 * ((t139 * t64 + (t46 + t260) * qJD(1)) * t139 + (t45 * qJD(1) + (-t95 * t208 - t97 * t210 + t224) * t137 + (-t63 + (t239 + t242) * qJD(4) - t180 * qJD(1)) * t139) * t137) + (t13 * t31 - t43 * t87 - t44 * t86) * t256 + t137 * t3 - t139 * t4 + (t137 * t46 - t45 * t139 + t15) * t212 + (t48 * t137 - t139 * t47 + t16) * t211; m(6) * (t19 * t51 + t20 * t50 + t22 * t53 + t23 * t52) + t49 + (-t21 + (t200 * t137 + t199 * t139) * qJD(4)) * t138 + ((t10 / 0.2e1 + t17 / 0.2e1) * t139 + (t9 / 0.2e1 + t18 / 0.2e1) * t137 + (-t199 * t137 + t200 * t139) * qJD(1)) * t136; m(6) * t14; m(6) * (t137 * t20 - t139 * t19 + (t137 * t53 + t139 * t52) * qJD(1)); m(6) * (t13 * t40 + t14 * t31 - t19 * t86 - t20 * t87 + t43 * t52 + t44 * t53) + (t12 * t244 + t16 * t194 - t2 / 0.2e1 + (qJD(1) * t30 - t9) * t252) * t139 + (t1 / 0.2e1 + t11 * t244 + t15 * t194 + (qJD(1) * t29 + t10) * t252) * t137 + (t3 * t250 + t4 * t253 + qJD(4) * t177 / 0.2e1 + (-t137 * t16 / 0.2e1 + t15 * t250) * qJD(1)) * t136; (t14 * t40 + t19 * t53 + t20 * t52) * t256 + (t21 * t138 - t49 + (t137 * t11 + t139 * t12 - t176 * t138) * qJD(4)) * t138 + (t139 * t1 + t137 * t2 - t138 * (t10 * t139 + t137 * t9) + (t176 * t136 - t54 * t138) * qJD(4) + (t139 * t11 - t137 * t12 + t138 * t177) * qJD(1)) * t136;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t57(1), t57(2), t57(4), t57(7), t57(11); t57(2), t57(3), t57(5), t57(8), t57(12); t57(4), t57(5), t57(6), t57(9), t57(13); t57(7), t57(8), t57(9), t57(10), t57(14); t57(11), t57(12), t57(13), t57(14), t57(15);];
Mq = res;
