% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPRPR5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR5_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_inertiaDJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:35:12
% EndTime: 2019-03-09 10:35:20
% DurationCPUTime: 2.88s
% Computational Cost: add. (5847->329), mult. (16671->623), div. (0->0), fcn. (16658->12), ass. (0->164)
t141 = sin(pkin(6));
t135 = t141 ^ 2;
t147 = sin(qJ(2));
t191 = qJD(2) * t147;
t213 = t135 * t191;
t199 = t141 * t147;
t144 = cos(pkin(6));
t210 = pkin(1) * t144;
t185 = t147 * t210;
t150 = cos(qJ(2));
t198 = t141 * t150;
t208 = pkin(8) + qJ(3);
t94 = t208 * t198 + t185;
t212 = -t94 * qJD(2) - qJD(3) * t199;
t142 = cos(pkin(12));
t132 = -t142 * pkin(5) - pkin(4);
t211 = 0.2e1 * t132;
t146 = sin(qJ(4));
t209 = pkin(4) * t146;
t207 = pkin(10) + qJ(5);
t139 = sin(pkin(12));
t140 = sin(pkin(11));
t143 = cos(pkin(11));
t100 = t140 * t199 - t143 * t198;
t133 = qJD(4) * t146;
t149 = cos(qJ(4));
t189 = qJD(4) * t149;
t190 = qJD(2) * t150;
t126 = t190 * t210;
t172 = t208 * t147;
t80 = t126 + (-qJD(2) * t172 + qJD(3) * t150) * t141;
t45 = t140 * t212 + t143 * t80;
t88 = (pkin(1) * t150 + pkin(2)) * t144 - t141 * t172;
t60 = t140 * t88 + t143 * t94;
t52 = t144 * pkin(9) + t60;
t177 = t141 * t191;
t125 = pkin(2) * t177;
t101 = (t140 * t150 + t143 * t147) * t141;
t95 = qJD(2) * t101;
t176 = t141 * t190;
t96 = -t140 * t177 + t143 * t176;
t64 = t95 * pkin(3) - t96 * pkin(9) + t125;
t182 = -pkin(2) * t150 - pkin(1);
t65 = t100 * pkin(3) - t101 * pkin(9) + t182 * t141;
t20 = t52 * t133 - t146 * t64 - t149 * t45 - t65 * t189;
t14 = t95 * qJ(5) + t100 * qJD(5) - t20;
t44 = t140 * t80 - t143 * t212;
t85 = t101 * t149 + t144 * t146;
t57 = qJD(4) * t85 + t96 * t146;
t84 = t101 * t146 - t144 * t149;
t58 = -qJD(4) * t84 + t96 * t149;
t24 = t57 * pkin(4) - t58 * qJ(5) - t85 * qJD(5) + t44;
t8 = t139 * t24 + t142 * t14;
t205 = t146 * t65 + t149 * t52;
t29 = t100 * qJ(5) + t205;
t59 = -t140 * t94 + t143 * t88;
t51 = -t144 * pkin(3) - t59;
t34 = t84 * pkin(4) - t85 * qJ(5) + t51;
t19 = t139 * t34 + t142 * t29;
t145 = sin(qJ(6));
t148 = cos(qJ(6));
t116 = t148 * t139 + t145 * t142;
t102 = t116 * t146;
t186 = qJD(6) * t148;
t187 = qJD(6) * t145;
t197 = t142 * t146;
t201 = t139 * t146;
t74 = t116 * t189 + t186 * t197 - t187 * t201;
t206 = -t102 * t57 - t74 * t84;
t21 = -t65 * t133 - t146 * t45 + t149 * t64 - t52 * t189;
t17 = -t95 * pkin(4) - t21;
t204 = t17 * t139;
t203 = t17 * t142;
t109 = -t146 * qJD(5) + (-qJ(5) * t149 + t209) * qJD(4);
t130 = t140 * pkin(2) + pkin(9);
t180 = t139 * t133;
t78 = t142 * t109 + t130 * t180;
t202 = t130 * t149;
t200 = t139 * t149;
t196 = t142 * t149;
t195 = t146 * t130;
t194 = t148 * t142;
t131 = -t143 * pkin(2) - pkin(3);
t162 = -t149 * pkin(4) - t146 * qJ(5);
t112 = t131 + t162;
t117 = t130 * t196;
t82 = t139 * t112 + t117;
t193 = t139 ^ 2 + t142 ^ 2;
t192 = t146 ^ 2 - t149 ^ 2;
t188 = qJD(5) * t149;
t184 = 0.2e1 * qJD(4) * t131;
t183 = t130 * t200;
t181 = t135 * t190;
t179 = t142 * t133;
t118 = t130 * t189;
t178 = t139 * t189;
t174 = t142 * t189;
t173 = t146 * t189;
t7 = -t139 * t14 + t142 * t24;
t18 = -t139 * t29 + t142 * t34;
t62 = -t100 * t142 + t85 * t139;
t63 = t100 * t139 + t85 * t142;
t35 = t145 * t63 + t148 * t62;
t40 = t58 * t139 - t95 * t142;
t41 = t95 * t139 + t58 * t142;
t15 = -qJD(6) * t35 - t145 * t40 + t148 * t41;
t36 = -t145 * t62 + t148 * t63;
t171 = t36 * t133 - t15 * t149;
t170 = -t146 * t52 + t149 * t65;
t169 = t193 * t149;
t168 = t193 * qJD(5);
t167 = 0.2e1 * t173;
t166 = t40 * pkin(10) - t8;
t165 = 0.2e1 * t168;
t164 = -t7 * t139 + t8 * t142;
t115 = t145 * t139 - t194;
t103 = t115 * t146;
t111 = t116 * qJD(6);
t73 = t111 * t146 + t145 * t178 - t148 * t174;
t163 = t103 * t57 + t73 * t84;
t11 = t84 * pkin(5) - t63 * pkin(10) + t18;
t12 = -t62 * pkin(10) + t19;
t4 = t145 * t11 + t148 * t12;
t161 = -t139 * t18 + t142 * t19;
t97 = t139 * t109;
t79 = -t130 * t179 + t97;
t160 = -t78 * t139 + t79 * t142;
t105 = t142 * t112;
t81 = t105 - t183;
t159 = -t139 * t81 + t142 * t82;
t72 = -pkin(10) * t197 + t105 + (-t130 * t139 - pkin(5)) * t149;
t75 = -pkin(10) * t201 + t82;
t39 = t145 * t72 + t148 * t75;
t158 = -qJ(5) * t57 - qJD(5) * t84;
t122 = t207 * t139;
t123 = t207 * t142;
t87 = -t145 * t122 + t148 * t123;
t157 = pkin(5) * t146 - pkin(10) * t196;
t30 = -t100 * pkin(4) - t170;
t16 = qJD(6) * t36 + t145 * t41 + t148 * t40;
t156 = -t35 * t133 + t149 * t16;
t155 = t100 * t189 + t146 * t95;
t154 = -t149 * t111 + t115 * t133;
t153 = -pkin(10) * t200 - t142 * t195;
t152 = t57 * pkin(5) - t41 * pkin(10) + t7;
t110 = t115 * qJD(6);
t108 = (pkin(5) * t139 + t130) * t146;
t107 = (-pkin(8) * t198 - t185) * qJD(2);
t106 = pkin(8) * t177 - t126;
t99 = pkin(5) * t178 + t118;
t86 = -t148 * t122 - t145 * t123;
t76 = t110 * t149 + t116 * t133;
t71 = -qJD(5) * t116 - qJD(6) * t87;
t70 = t122 * t186 - qJD(5) * t194 + (qJD(5) * t139 + qJD(6) * t123) * t145;
t67 = -t100 * t133 + t149 * t95;
t38 = -t145 * t75 + t148 * t72;
t27 = -t145 * t97 + t148 * t78 - t39 * qJD(6) + (-t145 * t153 + t148 * t157) * qJD(4);
t26 = t75 * t187 - t145 * (qJD(4) * t157 + t78) - t72 * t186 - t148 * (qJD(4) * t153 + t97);
t25 = t62 * pkin(5) + t30;
t9 = t40 * pkin(5) + t17;
t3 = t148 * t11 - t145 * t12;
t2 = -qJD(6) * t4 + t145 * t166 + t148 * t152;
t1 = -t11 * t186 + t12 * t187 - t145 * t152 + t148 * t166;
t5 = [0, 0, 0, 0.2e1 * t147 * t181, 0.2e1 * (-t147 ^ 2 + t150 ^ 2) * t135 * qJD(2), 0.2e1 * t144 * t176, -0.2e1 * t144 * t177, 0, -0.2e1 * pkin(1) * t213 + 0.2e1 * t107 * t144, -0.2e1 * pkin(1) * t181 + 0.2e1 * t106 * t144, -0.2e1 * t45 * t100 + 0.2e1 * t44 * t101 - 0.2e1 * t59 * t96 - 0.2e1 * t60 * t95, 0.2e1 * pkin(2) * t182 * t213 - 0.2e1 * t59 * t44 + 0.2e1 * t60 * t45, 0.2e1 * t85 * t58, -0.2e1 * t85 * t57 - 0.2e1 * t58 * t84, 0.2e1 * t58 * t100 + 0.2e1 * t85 * t95, -0.2e1 * t57 * t100 - 0.2e1 * t84 * t95, 0.2e1 * t100 * t95, 0.2e1 * t21 * t100 + 0.2e1 * t170 * t95 + 0.2e1 * t44 * t84 + 0.2e1 * t51 * t57, 0.2e1 * t20 * t100 - 0.2e1 * t205 * t95 + 0.2e1 * t44 * t85 + 0.2e1 * t51 * t58, 0.2e1 * t17 * t62 + 0.2e1 * t18 * t57 + 0.2e1 * t30 * t40 + 0.2e1 * t7 * t84, 0.2e1 * t17 * t63 - 0.2e1 * t19 * t57 + 0.2e1 * t30 * t41 - 0.2e1 * t8 * t84, -0.2e1 * t18 * t41 - 0.2e1 * t19 * t40 - 0.2e1 * t8 * t62 - 0.2e1 * t7 * t63, 0.2e1 * t30 * t17 + 0.2e1 * t18 * t7 + 0.2e1 * t19 * t8, 0.2e1 * t36 * t15, -0.2e1 * t15 * t35 - 0.2e1 * t36 * t16, 0.2e1 * t15 * t84 + 0.2e1 * t36 * t57, -0.2e1 * t16 * t84 - 0.2e1 * t35 * t57, 0.2e1 * t84 * t57, 0.2e1 * t25 * t16 + 0.2e1 * t2 * t84 + 0.2e1 * t3 * t57 + 0.2e1 * t9 * t35, 0.2e1 * t1 * t84 + 0.2e1 * t25 * t15 + 0.2e1 * t9 * t36 - 0.2e1 * t4 * t57; 0, 0, 0, 0, 0, t176, -t177, 0, t107, t106 (-t140 * t95 - t143 * t96) * pkin(2) (t140 * t45 - t143 * t44) * pkin(2), t58 * t146 + t85 * t189, -t146 * t57 + t58 * t149 + (-t146 * t85 - t149 * t84) * qJD(4), t155, t67, 0, -t95 * t195 + t131 * t57 - t44 * t149 + (-t100 * t202 + t146 * t51) * qJD(4), -t95 * t202 + t131 * t58 + t44 * t146 + (t100 * t195 + t149 * t51) * qJD(4), -t7 * t149 + t81 * t57 + t78 * t84 + (t130 * t40 + t204) * t146 + (t146 * t18 + (t130 * t62 + t139 * t30) * t149) * qJD(4), t8 * t149 - t82 * t57 - t79 * t84 + (t130 * t41 + t203) * t146 + (-t146 * t19 + (t130 * t63 + t142 * t30) * t149) * qJD(4), -t82 * t40 - t81 * t41 - t79 * t62 - t78 * t63 + (-t139 * t8 - t142 * t7) * t146 + (-t139 * t19 - t142 * t18) * t189, t18 * t78 + t19 * t79 + t7 * t81 + t8 * t82 + (t146 * t17 + t30 * t189) * t130, -t15 * t103 - t36 * t73, -t15 * t102 + t103 * t16 + t73 * t35 - t36 * t74, -t163 + t171, t156 + t206, t133 * t84 - t57 * t149, t9 * t102 + t108 * t16 + t133 * t3 - t2 * t149 + t25 * t74 + t27 * t84 + t99 * t35 + t38 * t57, -t1 * t149 - t9 * t103 + t108 * t15 - t133 * t4 - t25 * t73 + t26 * t84 + t99 * t36 - t39 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t167, -0.2e1 * t192 * qJD(4), 0, 0, 0, t146 * t184, t149 * t184, -0.2e1 * t78 * t149 + 0.2e1 * (t81 + 0.2e1 * t183) * t133, 0.2e1 * t79 * t149 + 0.2e1 * (-t82 + 0.2e1 * t117) * t133, 0.2e1 * (-t139 * t79 - t142 * t78) * t146 + 0.2e1 * (-t139 * t82 - t142 * t81) * t189, 0.2e1 * t130 ^ 2 * t173 + 0.2e1 * t81 * t78 + 0.2e1 * t82 * t79, 0.2e1 * t103 * t73, 0.2e1 * t73 * t102 + 0.2e1 * t103 * t74, -0.2e1 * t103 * t133 + 0.2e1 * t73 * t149, -0.2e1 * t102 * t133 + 0.2e1 * t149 * t74, -0.2e1 * t173, 0.2e1 * t99 * t102 + 0.2e1 * t108 * t74 + 0.2e1 * t133 * t38 - 0.2e1 * t27 * t149, -0.2e1 * t99 * t103 - 0.2e1 * t108 * t73 - 0.2e1 * t133 * t39 - 0.2e1 * t26 * t149; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t125, 0, 0, 0, 0, 0, t67, -t155, -t57 * t201 - t149 * t40 + (t146 * t62 - t84 * t200) * qJD(4), -t57 * t197 - t149 * t41 + (t146 * t63 - t84 * t196) * qJD(4) (t139 * t41 - t142 * t40) * t146 + (t139 * t63 - t142 * t62) * t189, -t17 * t149 + t164 * t146 + (t146 * t30 + t149 * t161) * qJD(4), 0, 0, 0, 0, 0, -t156 + t206, t163 + t171; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t160 * t146 + (t192 * t130 + t159 * t149) * qJD(4), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-0.1e1 + t193) * t167, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, -t57, t95, t21, t20, -pkin(4) * t40 + t139 * t158 - t203, -pkin(4) * t41 + t142 * t158 + t204 (-qJ(5) * t40 - qJD(5) * t62 + t8) * t142 + (qJ(5) * t41 + qJD(5) * t63 - t7) * t139, -t17 * pkin(4) + qJ(5) * t164 + qJD(5) * t161, -t36 * t110 + t15 * t116, t110 * t35 - t36 * t111 - t15 * t115 - t116 * t16, -t110 * t84 + t116 * t57, -t111 * t84 - t115 * t57, 0, t25 * t111 + t9 * t115 + t132 * t16 + t86 * t57 + t71 * t84, -t25 * t110 + t9 * t116 + t132 * t15 - t87 * t57 + t70 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t189, -t133, 0, -t118, t130 * t133, t139 * t188 + (t139 * t162 - t117) * qJD(4), t142 * t188 + (t142 * t162 + t183) * qJD(4), t160, -pkin(4) * t118 + qJ(5) * t160 + qJD(5) * t159, t103 * t110 - t73 * t116, t110 * t102 + t103 * t111 + t73 * t115 - t116 * t74, t76, -t154, 0, t108 * t111 + t99 * t115 + t132 * t74 + t133 * t86 - t71 * t149, -t108 * t110 + t99 * t116 - t132 * t73 - t133 * t87 - t70 * t149; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t133, -t189, -t179, t180, qJD(4) * t169, t146 * t168 + (qJ(5) * t169 - t209) * qJD(4), 0, 0, 0, 0, 0, t154, t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t165, qJ(5) * t165, -0.2e1 * t116 * t110, 0.2e1 * t110 * t115 - 0.2e1 * t116 * t111, 0, 0, 0, t111 * t211, -t110 * t211; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t41, 0, t17, 0, 0, 0, 0, 0, t16, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t178, t174, 0, t118, 0, 0, 0, 0, 0, t74, -t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t133, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t111, -t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t16, t57, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73, -t74, t133, t27, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110, -t111, 0, t71, t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t5;
