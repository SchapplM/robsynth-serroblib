% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RRRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRRRPR3_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR3_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR3_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR3_inertiaDJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:04:00
% EndTime: 2019-03-09 22:04:09
% DurationCPUTime: 3.69s
% Computational Cost: add. (7904->282), mult. (17075->452), div. (0->0), fcn. (16748->8), ass. (0->156)
t195 = sin(qJ(3));
t203 = -pkin(8) - pkin(7);
t209 = t195 * t203;
t213 = -t195 * pkin(9) + t209;
t197 = cos(qJ(3));
t208 = t197 * t203;
t214 = -t197 * pkin(9) + t208;
t88 = sin(qJ(2));
t90 = cos(qJ(2));
t102 = t213 * t90 + t214 * t88;
t103 = t213 * t88 - t214 * t90;
t196 = cos(qJ(4));
t87 = sin(qJ(4));
t30 = -t196 * t102 + t87 * t103;
t111 = -t88 * t208 - t209 * t90;
t212 = t111 * qJD(3);
t86 = sin(qJ(6));
t84 = t86 ^ 2;
t89 = cos(qJ(6));
t85 = t89 ^ 2;
t182 = t84 - t85;
t206 = t182 * qJD(6);
t211 = t208 * t90 - t88 * t209;
t149 = t196 * t195;
t200 = t87 * pkin(3);
t80 = qJD(4) * t200;
t44 = t80 + (qJD(3) + qJD(4)) * pkin(2) * (t197 * t87 + t149);
t37 = t80 + t44;
t126 = t195 * t90 + t197 * t88;
t116 = t126 * qJD(3);
t104 = -qJD(2) * t126 - t116;
t159 = t195 * qJD(3);
t160 = t197 * qJD(3);
t161 = t195 * t88;
t205 = -(t197 * qJD(2) + t160) * t90 + qJD(2) * t161 + t88 * t159;
t125 = -t197 * t90 + t161;
t50 = -t87 * t125 + t196 * t126;
t25 = t50 * qJD(4) - t196 * t104 - t205 * t87;
t175 = t88 * qJD(2);
t169 = pkin(2) * t175;
t49 = t196 * t125 + t87 * t126;
t24 = t49 * qJD(4) - t87 * t104 + t196 * t205;
t9 = -pkin(3) * t104 + t25 * pkin(4) + t24 * qJ(5) - t50 * qJD(5) + t169;
t98 = t50 * pkin(5) + t30;
t210 = -t25 * pkin(10) - qJD(6) * t98 - t9;
t31 = t87 * t102 + t196 * t103;
t207 = t211 * qJD(3);
t168 = t197 * pkin(2);
t150 = t168 + pkin(3);
t128 = t196 * t150;
t157 = pkin(2) * t160;
t43 = (t195 * qJD(4) + t159) * pkin(2) * t87 - qJD(4) * t128 - t196 * t157;
t92 = 0.2e1 * qJD(5);
t204 = pkin(4) + pkin(10);
t178 = qJD(6) * t86;
t79 = -t90 * pkin(2) - pkin(1);
t54 = pkin(3) * t125 + t79;
t105 = -t50 * qJ(5) + t54;
t22 = t204 * t49 + t105;
t136 = qJD(2) * t209;
t137 = qJD(2) * t208;
t94 = pkin(9) * t104 + t136 * t90 + t137 * t88 - t212;
t95 = pkin(9) * t205 - t136 * t88 + t137 * t90 + t207;
t11 = t31 * qJD(4) - t196 * t95 + t87 * t94;
t93 = -t24 * pkin(5) + t11;
t2 = t22 * t178 + t210 * t89 - t86 * t93;
t201 = t2 * t86;
t10 = t30 * qJD(4) - t196 * t94 - t87 * t95;
t8 = -t25 * pkin(5) - t10;
t6 = t8 * t86;
t7 = t8 * t89;
t199 = t88 * pkin(2);
t177 = qJD(6) * t89;
t23 = -t49 * pkin(5) + t31;
t198 = t23 * t177 + t6;
t194 = t25 * t86;
t193 = t25 * t89;
t192 = t30 * t44;
t41 = -qJD(5) + t43;
t61 = pkin(2) * t149 + t87 * t150;
t56 = qJ(5) + t61;
t191 = t41 * t56;
t190 = t44 * t50;
t189 = t49 * t25;
t166 = t196 * pkin(3);
t155 = qJD(4) * t166;
t71 = t155 + qJD(5);
t76 = qJ(5) + t200;
t188 = t71 * t76;
t187 = t56 * t177 - t41 * t86;
t186 = t76 * t177 + t71 * t86;
t173 = qJ(5) * qJD(6);
t183 = qJD(5) * t86 + t89 * t173;
t181 = t84 + t85;
t180 = qJD(6) * t23;
t176 = qJD(6) * t204;
t174 = t90 * qJD(2);
t172 = 0.2e1 * t189;
t19 = -0.2e1 * t50 * t24;
t171 = -0.2e1 * pkin(1) * qJD(2);
t170 = t86 * t193;
t165 = t195 * pkin(2);
t163 = t88 * t174;
t162 = t86 * t177;
t32 = t181 * t44;
t48 = t49 ^ 2;
t158 = t48 * t162;
t156 = pkin(2) * t159;
t78 = -t166 - pkin(4);
t148 = -t10 * t31 + t11 * t30;
t14 = -t86 * t22 + t89 * t98;
t15 = t89 * t22 + t86 * t98;
t147 = t14 * t89 + t15 * t86;
t146 = -t14 * t86 + t15 * t89;
t145 = -t41 * t76 + t56 * t71;
t59 = t181 * t80;
t140 = -qJ(5) * t25 - qJD(5) * t49;
t139 = -qJ(5) * t41 + qJD(5) * t56;
t138 = qJ(5) * t71 + qJD(5) * t76;
t135 = t50 * t177 - t24 * t86;
t17 = -t50 * t178 - t24 * t89;
t134 = t49 * t177 + t194;
t133 = t49 * t178 - t193;
t60 = -t87 * t165 + t128;
t57 = -pkin(4) - t60;
t55 = -pkin(10) + t57;
t132 = qJD(6) * (t49 * t56 - t50 * t55);
t75 = -pkin(10) + t78;
t131 = qJD(6) * (t49 * t76 - t50 * t75);
t130 = qJD(6) * (qJ(5) * t49 + t204 * t50);
t129 = 0.2e1 * t24 * t49 - 0.2e1 * t25 * t50;
t127 = -t25 * t56 + t41 * t49 + t190;
t124 = t204 * t24 + t140;
t120 = t79 * t126;
t117 = -t25 * t76 - t49 * t71 + t50 * t80;
t115 = -t24 * t55 + t127;
t112 = -t24 * t75 + t117;
t3 = -t22 * t177 + t210 * t86 + t89 * t93;
t110 = qJD(6) * t146 + t3 * t89 - t201;
t109 = 0.2e1 * t10 * t49 + 0.2e1 * t11 * t50 - 0.2e1 * t24 * t30 - 0.2e1 * t25 * t31;
t106 = -t43 + t92;
t83 = qJ(5) * t92;
t82 = qJD(5) * t89;
t70 = -0.2e1 * t162;
t69 = 0.2e1 * t162;
t65 = t71 * t89;
t63 = 0.2e1 * t206;
t40 = pkin(3) * t116 + (pkin(3) * t126 + t199) * qJD(2);
t39 = t41 * t89;
t34 = qJD(2) * t211 + t207;
t33 = qJD(2) * t111 + t212;
t29 = t49 * pkin(4) + t105;
t16 = -t206 * t49 + t170;
t13 = -0.4e1 * t49 * t162 - t182 * t25;
t1 = -t14 * t178 - t201 + (qJD(6) * t15 + t3) * t89;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t163, 0.2e1 * (-t88 ^ 2 + t90 ^ 2) * qJD(2), 0, -0.2e1 * t163, 0, 0, t88 * t171, t90 * t171, 0, 0, -0.2e1 * t126 * t205, 0.2e1 * t104 * t126 + 0.2e1 * t125 * t205, 0, -0.2e1 * t125 * t104, 0, 0, 0.2e1 * qJD(3) * t120 + 0.2e1 * (t125 * t199 + t120) * qJD(2), 0.2e1 * t126 * t169 - 0.2e1 * t205 * t79, -0.2e1 * t104 * t211 - 0.2e1 * t111 * t205 + 0.2e1 * t125 * t33 - 0.2e1 * t126 * t34, -0.2e1 * t111 * t34 + 0.2e1 * t79 * t169 + 0.2e1 * t211 * t33, t19, t129, 0, t172, 0, 0, 0.2e1 * t25 * t54 + 0.2e1 * t40 * t49, -0.2e1 * t24 * t54 + 0.2e1 * t40 * t50, t109, 0.2e1 * t40 * t54 + 0.2e1 * t148, 0, 0, 0, t19, t129, t172, t109, -0.2e1 * t25 * t29 - 0.2e1 * t49 * t9, 0.2e1 * t24 * t29 - 0.2e1 * t50 * t9, 0.2e1 * t29 * t9 + 0.2e1 * t148, 0.2e1 * t84 * t189 + 0.2e1 * t158, 0.4e1 * t49 * t170 - 0.2e1 * t48 * t206, 0.2e1 * t135 * t49 + 0.2e1 * t50 * t194, 0.2e1 * t85 * t189 - 0.2e1 * t158, 0.2e1 * t17 * t49 + 0.2e1 * t50 * t193, t19, 0.2e1 * t133 * t23 - 0.2e1 * t14 * t24 + 0.2e1 * t3 * t50 - 0.2e1 * t49 * t7, 0.2e1 * t134 * t23 + 0.2e1 * t15 * t24 + 0.2e1 * t2 * t50 + 0.2e1 * t49 * t6, 0.2e1 * t146 * t25 + 0.2e1 * (-qJD(6) * t147 - t2 * t89 - t3 * t86) * t49, 0.2e1 * t14 * t3 - 0.2e1 * t15 * t2 + 0.2e1 * t23 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t174, 0, -t175, 0, -pkin(7) * t174, pkin(7) * t175, 0, 0, 0, 0, -t205, 0, t104, 0, t34, t33, t205 * t168 - t125 * t157 + (t104 + t116) * t165 (t34 * t197 - t33 * t195 + (t111 * t195 - t197 * t211) * qJD(3)) * pkin(2), 0, 0, -t24, 0, -t25, 0, -t11, t10, t24 * t60 - t25 * t61 + t43 * t49 + t190, -t10 * t61 - t11 * t60 - t31 * t43 + t192, 0, t24, t25, 0, 0, 0, -t24 * t57 + t127, t11, -t10, -t10 * t56 + t11 * t57 - t31 * t41 + t192, t16, t13, t17, -t16, -t135, 0, t115 * t89 + t132 * t86 + t198, t7 + t89 * t132 + (-t115 - t180) * t86, -t1, t110 * t55 + t147 * t44 - t23 * t41 + t56 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t156, -0.2e1 * t157, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t44, 0.2e1 * t43, 0, -0.2e1 * t43 * t61 - 0.2e1 * t44 * t60, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t44, -0.2e1 * t41, 0.2e1 * t44 * t57 - 0.2e1 * t191, t70, t63, 0, t69, 0, 0, 0.2e1 * t187, -0.2e1 * t56 * t178 - 0.2e1 * t39, -0.2e1 * t32, 0.2e1 * t32 * t55 - 0.2e1 * t191; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t205, 0, t104, 0, t34, t33, 0, 0, 0, 0, -t24, 0, -t25, 0, -t11, t10 (t196 * t24 - t25 * t87 + (-t196 * t49 + t50 * t87) * qJD(4)) * pkin(3) (-t196 * t11 - t10 * t87 + (t196 * t31 + t30 * t87) * qJD(4)) * pkin(3), 0, t24, t25, 0, 0, 0, -t24 * t78 + t117, t11, -t10, -t10 * t76 + t11 * t78 + t30 * t80 + t31 * t71, t16, t13, t17, -t16, -t135, 0, t112 * t89 + t131 * t86 + t198, t7 + t89 * t131 + (-t112 - t180) * t86, -t1, t110 * t75 + t147 * t80 + t23 * t71 + t76 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t156, -t157, 0, 0, 0, 0, 0, 0, 0, 0, -t37, -t155 + t43, 0 (-t196 * t44 - t43 * t87 + (t196 * t61 - t60 * t87) * qJD(4)) * pkin(3), 0, 0, 0, 0, 0, 0, 0, t37, t155 + t106, t44 * t78 + t57 * t80 + t145, t70, t63, 0, t69, 0, 0, t186 + t187, -t39 + t65 + (-t56 - t76) * t178, -t37 * t181, t32 * t75 + t55 * t59 + t145; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t80, -0.2e1 * t155, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t80, 0.2e1 * t71, 0.2e1 * t78 * t80 + 0.2e1 * t188, t70, t63, 0, t69, 0, 0, 0.2e1 * t186, -0.2e1 * t76 * t178 + 0.2e1 * t65, -0.2e1 * t59, 0.2e1 * t59 * t75 + 0.2e1 * t188; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, 0, -t25, 0, -t11, t10, 0, 0, 0, t24, t25, 0, 0, 0, pkin(4) * t24 + t140, t11, -t10, -pkin(4) * t11 - qJ(5) * t10 + qJD(5) * t31, t16, t13, t17, -t16, -t135, 0, t124 * t89 + t130 * t86 + t198, t7 + t89 * t130 + (-t124 - t180) * t86, -t1, qJ(5) * t8 + qJD(5) * t23 - t110 * t204; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, t43, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t106, -pkin(4) * t44 + t139, t70, t63, 0, t69, 0, 0, t183 + t187, -t39 + t82 + (-qJ(5) - t56) * t178, -t32, -t204 * t32 + t139; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80, -t155, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, t92 + t155, -pkin(4) * t80 + t138, t70, t63, 0, t69, 0, 0, t183 + t186, t65 + t82 + (-qJ(5) - t76) * t178, -t59, -t204 * t59 + t138; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, t83, t70, t63, 0, t69, 0, 0, 0.2e1 * t183, -0.2e1 * t86 * t173 + 0.2e1 * t82, 0, t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, 0, 0, t11, 0, 0, 0, 0, 0, 0, t17, -t135, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134, 0, -t133, -t24, t3, t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t178, 0, -t177, 0, -t55 * t178 + t44 * t89, -t55 * t177 - t44 * t86, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t178, 0, -t177, 0, -t75 * t178 + t89 * t80, -t75 * t177 - t86 * t80, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t178, 0, -t177, 0, t86 * t176, t89 * t176, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t178, -t177, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t4;