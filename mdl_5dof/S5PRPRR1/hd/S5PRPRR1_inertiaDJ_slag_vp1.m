% Calculate time derivative of joint inertia matrix for
% S5PRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR1_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR1_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR1_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR1_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR1_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR1_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:42:38
% EndTime: 2019-12-05 15:42:46
% DurationCPUTime: 2.74s
% Computational Cost: add. (6353->289), mult. (4442->425), div. (0->0), fcn. (3422->8), ass. (0->172)
t119 = pkin(8) + qJ(2);
t112 = sin(t119);
t114 = cos(t119);
t118 = pkin(9) + qJ(4);
t115 = qJ(5) + t118;
t106 = sin(t115);
t204 = rSges(6,2) * t106;
t107 = cos(t115);
t206 = rSges(6,1) * t107;
t155 = -t204 + t206;
t68 = t112 * rSges(6,3) + t155 * t114;
t111 = sin(t118);
t113 = cos(t118);
t94 = rSges(5,1) * t111 + rSges(5,2) * t113;
t132 = qJD(4) * t94;
t223 = t112 * t132;
t120 = qJD(4) + qJD(5);
t203 = rSges(6,2) * t107;
t89 = rSges(6,1) * t106 + t203;
t133 = t89 * t120;
t170 = qJD(4) * t111;
t167 = pkin(4) * t170;
t221 = t167 + t133;
t138 = Icges(6,5) * t107 - Icges(6,6) * t106;
t61 = -Icges(6,3) * t114 + t138 * t112;
t220 = qJD(2) * t61;
t189 = Icges(5,4) * t113;
t142 = -Icges(5,2) * t111 + t189;
t72 = Icges(5,6) * t112 + t142 * t114;
t190 = Icges(5,4) * t111;
t145 = Icges(5,1) * t113 - t190;
t74 = Icges(5,5) * t112 + t145 * t114;
t148 = t111 * t72 - t113 * t74;
t219 = t148 * t112;
t71 = -Icges(5,6) * t114 + t142 * t112;
t73 = -Icges(5,5) * t114 + t145 * t112;
t149 = t111 * t71 - t113 * t73;
t218 = t149 * t114;
t187 = Icges(6,4) * t107;
t140 = -Icges(6,2) * t106 + t187;
t64 = Icges(6,6) * t112 + t140 * t114;
t188 = Icges(6,4) * t106;
t143 = Icges(6,1) * t107 - t188;
t66 = Icges(6,5) * t112 + t143 * t114;
t153 = t106 * t64 - t107 * t66;
t217 = t153 * t112;
t63 = -Icges(6,6) * t114 + t140 * t112;
t65 = -Icges(6,5) * t114 + t143 * t112;
t154 = t106 * t63 - t107 * t65;
t216 = t154 * t114;
t123 = -pkin(6) - qJ(3);
t117 = -pkin(7) + t123;
t173 = t117 - t123;
t215 = t173 * t114;
t122 = cos(pkin(9));
t108 = t122 * pkin(3) + pkin(2);
t205 = rSges(5,2) * t111;
t207 = rSges(5,1) * t113;
t156 = -t205 + t207;
t135 = -t108 - t156;
t50 = (rSges(5,3) - t123) * t114 + t135 * t112;
t105 = t112 * rSges(5,3);
t191 = t114 * t207 + t105;
t51 = -t112 * t123 + (t108 - t205) * t114 + t191;
t214 = t112 * t51 + t114 * t50;
t87 = Icges(6,2) * t107 + t188;
t88 = Icges(6,1) * t106 + t187;
t152 = t106 * t87 - t107 * t88;
t213 = t152 * qJD(2) + t138 * t120;
t139 = Icges(5,5) * t113 - Icges(5,6) * t111;
t69 = -Icges(5,3) * t114 + t139 * t112;
t137 = rSges(4,1) * t122 - rSges(4,2) * sin(pkin(9)) + pkin(2);
t194 = rSges(4,3) + qJ(3);
t55 = t194 * t112 + t137 * t114;
t212 = 2 * m(5);
t211 = 2 * m(6);
t109 = t112 ^ 2;
t110 = t114 ^ 2;
t210 = m(5) * t94;
t209 = t112 / 0.2e1;
t208 = -t114 / 0.2e1;
t67 = -rSges(6,3) * t114 + t155 * t112;
t31 = t112 * t67 + t114 * t68;
t202 = t111 * t73;
t201 = t111 * t74;
t199 = t113 * t71;
t198 = t113 * t72;
t196 = t120 * t87;
t195 = t120 * t88;
t193 = rSges(6,3) - t117;
t171 = qJD(2) * t114;
t172 = qJD(2) * t112;
t192 = rSges(6,3) * t171 + t172 * t204;
t62 = Icges(6,3) * t112 + t138 * t114;
t180 = qJD(2) * t62;
t70 = Icges(5,3) * t112 + t139 * t114;
t179 = qJD(2) * t70;
t178 = t106 * t120;
t177 = t107 * t120;
t176 = t112 * t117;
t175 = t114 * t120;
t174 = t109 + t110;
t169 = qJD(4) * t113;
t168 = t112 * (t68 * qJD(2) - t112 * t133) + t114 * (-t175 * t203 + (-t106 * t175 - t107 * t172) * rSges(6,1) + t192) + t67 * t171;
t166 = t114 * t205;
t165 = t111 * t172;
t164 = -pkin(4) * t111 - t89;
t36 = -t63 * qJD(2) - t114 * t196;
t161 = t120 * t66 + t36;
t37 = t64 * qJD(2) - t112 * t196;
t160 = t120 * t65 + t37;
t38 = -t65 * qJD(2) - t114 * t195;
t159 = -t120 * t64 + t38;
t39 = t66 * qJD(2) - t112 * t195;
t158 = t120 * t63 - t39;
t13 = -t154 * t112 - t61 * t114;
t14 = -t114 * t62 - t217;
t15 = t112 * t61 - t216;
t16 = t62 * t112 - t153 * t114;
t86 = Icges(6,5) * t106 + Icges(6,6) * t107;
t129 = t86 * t120;
t34 = -t114 * t129 - t220;
t35 = -t112 * t129 + t180;
t157 = -t114 * ((t114 * t35 + (t14 + t216) * qJD(2)) * t114 + (t13 * qJD(2) + (-t106 * t36 + t107 * t38 - t64 * t177 - t66 * t178 + t180) * t112 + (-t34 + t158 * t107 + t160 * t106 + (-t153 - t61) * qJD(2)) * t114) * t112) + t112 * ((t112 * t34 + (t15 + t217) * qJD(2)) * t112 + (t16 * qJD(2) + (t106 * t37 - t107 * t39 + t63 * t177 + t65 * t178 - t220) * t114 + (-t35 + t159 * t107 - t161 * t106 + (-t154 + t62) * qJD(2)) * t112) * t114) + (t112 * t14 - t114 * t13) * t172 + (t112 * t16 - t114 * t15) * t171;
t97 = pkin(4) * t113 + t108;
t136 = -t155 - t97;
t46 = t136 * t112 + t193 * t114;
t81 = t114 * t97;
t47 = t68 + t81 - t176;
t147 = t112 * t47 + t114 * t46;
t59 = t164 * t112;
t60 = t164 * t114;
t146 = t112 * t59 + t114 * t60;
t144 = Icges(5,1) * t111 + t189;
t141 = Icges(5,2) * t113 + t190;
t78 = t140 * t120;
t79 = t143 * t120;
t124 = qJD(2) * t86 + (t79 - t196) * t107 + (-t78 - t195) * t106;
t134 = (t159 * t106 + t161 * t107 + t213 * t112 + t124 * t114) * t209 + (-t158 * t106 + t160 * t107 + t124 * t112 - t213 * t114) * t208 + (t106 * t65 + t107 * t63 - t152 * t112 - t114 * t86) * t172 / 0.2e1 + (t106 * t66 + t107 * t64 + t112 * t86 - t152 * t114) * t171 / 0.2e1;
t128 = qJD(4) * t144;
t127 = qJD(4) * t141;
t126 = qJD(4) * (-Icges(5,5) * t111 - Icges(5,6) * t113);
t125 = rSges(5,2) * t165 + rSges(5,3) * t171 - t114 * t132;
t54 = -t137 * t112 + t194 * t114;
t103 = qJD(3) * t114;
t102 = qJD(3) * t112;
t99 = t123 * t172;
t85 = t156 * qJD(4);
t80 = t155 * t120;
t76 = -t166 + t191;
t75 = -t114 * rSges(5,3) + t156 * t112;
t53 = -t114 * t108 - t173 * t112 + t81;
t52 = t215 + t112 * (-t108 + t97);
t49 = -t55 * qJD(2) + t103;
t48 = t54 * qJD(2) + t102;
t41 = t112 * t126 + t179;
t40 = -qJD(2) * t69 + t114 * t126;
t30 = -t89 * t171 - t112 * t80 + (-t111 * t171 - t112 * t169) * pkin(4);
t29 = t89 * t172 - t114 * t80 + (-t114 * t169 + t165) * pkin(4);
t24 = t103 + t99 + t223 + (t135 * t114 - t105) * qJD(2);
t23 = t102 + (-t114 * t123 + (-t108 - t207) * t112) * qJD(2) + t125;
t22 = t103 + t221 * t112 + (-t193 * t112 + t136 * t114) * qJD(2);
t21 = t102 + (-t97 - t206) * t172 + (-qJD(2) * t117 - t221) * t114 + t192;
t20 = t112 * t70 - t148 * t114;
t19 = t112 * t69 - t218;
t18 = -t114 * t70 - t219;
t17 = -t149 * t112 - t114 * t69;
t12 = t112 * t52 + t114 * t53 + t31;
t11 = (qJD(2) * t75 + t125) * t114 + (-t223 + (-t166 - t76 + t105) * qJD(2)) * t112;
t10 = -t68 * t172 + t168;
t3 = t112 * (-t112 * t167 + t99) - t110 * t167 + ((t52 - t215) * t114 + (-t53 - t68 - t176) * t112) * qJD(2) + t168;
t1 = [0; 0; t88 * t177 + t106 * t79 - t87 * t178 + t107 * t78 + (t21 * t47 + t22 * t46) * t211 + (t23 * t51 + t24 * t50) * t212 + 0.2e1 * m(4) * (t48 * t55 + t49 * t54) + (t145 - t141) * t170 + (t144 + t142) * t169; 0; m(6) * (t147 * qJD(2) + t112 * t22 - t114 * t21) + m(5) * (qJD(2) * t214 + t112 * t24 - t114 * t23) + m(4) * (t112 * t49 - t114 * t48 + (t112 * t55 + t114 * t54) * qJD(2)); 0; m(5) * t11 + m(6) * t3; (-t148 * qJD(4) + t111 * (-t73 * qJD(2) - t114 * t128) + t113 * (-t71 * qJD(2) - t114 * t127)) * t209 + (-t149 * qJD(4) + t111 * (t74 * qJD(2) - t112 * t128) + t113 * (t72 * qJD(2) - t112 * t127)) * t208 + m(6) * (t59 * t21 + t22 * t60 + t29 * t46 + t30 * t47) + m(5) * ((-t112 * t23 - t114 * t24) * t94 - t214 * t85) + (t110 / 0.2e1 + t109 / 0.2e1) * t139 * qJD(4) + ((t201 / 0.2e1 + t198 / 0.2e1 - t51 * t210) * t114 + (t202 / 0.2e1 + t199 / 0.2e1 + t50 * t210) * t112) * qJD(2) + t134; m(6) * (t146 * qJD(2) + t112 * t29 - t114 * t30); ((t112 * t75 + t114 * t76) * t11 + t174 * t94 * t85) * t212 + (t112 * t20 - t114 * t19) * t171 + t112 * ((t112 * t40 + (t19 + t219) * qJD(2)) * t112 + (t20 * qJD(2) + (t71 * t169 + t73 * t170) * t114 + (-t41 + (-t198 - t201) * qJD(4) + (-t149 + t70) * qJD(2)) * t112) * t114) + (t112 * t18 - t114 * t17) * t172 - t114 * ((t114 * t41 + (t18 + t218) * qJD(2)) * t114 + (t17 * qJD(2) + (-t72 * t169 - t74 * t170 + t179) * t112 + (-t40 + (t199 + t202) * qJD(4) - t148 * qJD(2)) * t114) * t112) + (t12 * t3 + t29 * t60 + t59 * t30) * t211 + t157; m(6) * t10; m(6) * (-t147 * t80 + (-t112 * t21 - t114 * t22 + (t112 * t46 - t114 * t47) * qJD(2)) * t89) + t134; 0; m(6) * (t10 * t12 + t31 * t3 - t146 * t80 + (-t112 * t30 - t114 * t29 + (t112 * t60 - t114 * t59) * qJD(2)) * t89) + t157; (t174 * t89 * t80 + t31 * t10) * t211 + t157;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
