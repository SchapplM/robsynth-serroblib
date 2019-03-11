% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RPRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRRPR7_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR7_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR7_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_inertiaDJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:21:01
% EndTime: 2019-03-09 05:21:08
% DurationCPUTime: 2.53s
% Computational Cost: add. (5212->201), mult. (10152->341), div. (0->0), fcn. (9816->8), ass. (0->128)
t142 = cos(pkin(10));
t87 = sin(qJ(4));
t88 = sin(qJ(3));
t90 = cos(qJ(4));
t91 = cos(qJ(3));
t66 = t87 * t91 + t90 * t88;
t67 = -t87 * t88 + t90 * t91;
t85 = sin(pkin(10));
t100 = t142 * t66 + t85 * t67;
t164 = qJD(3) + qJD(4);
t44 = t164 * t66;
t45 = t164 * t67;
t101 = t142 * t45 - t85 * t44;
t166 = t100 * t101;
t27 = t142 * t44 + t45 * t85;
t40 = -t142 * t67 + t66 * t85;
t176 = t40 * t27;
t99 = 0.2e1 * t176 + 0.2e1 * t166;
t182 = (-t101 * t85 + t142 * t27) * pkin(4);
t124 = t142 * t87;
t143 = pkin(3) * qJD(4);
t54 = (t85 * t90 + t124) * t143;
t153 = t40 * t54;
t133 = t90 * t143;
t134 = t87 * t143;
t55 = t142 * t133 - t85 * t134;
t112 = -t100 * t55 - t153;
t80 = pkin(3) * t90 + pkin(4);
t56 = -t85 * t87 * pkin(3) + t142 * t80;
t57 = pkin(3) * t124 + t85 * t80;
t181 = -t101 * t57 + t27 * t56 + t112;
t37 = t40 ^ 2;
t179 = (-t44 * t90 + t45 * t87 + (t66 * t90 - t67 * t87) * qJD(4)) * pkin(3);
t52 = -pkin(5) - t56;
t178 = t27 * t52;
t75 = -t142 * pkin(4) - pkin(5);
t177 = t27 * t75;
t158 = -pkin(1) - pkin(7);
t136 = pkin(8) - t158;
t120 = qJD(4) * t136;
t122 = t91 * t136;
t175 = qJD(3) * t122 + t91 * t120;
t139 = t88 * qJD(3);
t174 = t88 * t120 + t136 * t139;
t123 = t88 * t136;
t47 = -t87 * t122 - t90 * t123;
t33 = -qJ(5) * t66 + t47;
t46 = -t90 * t122 + t87 * t123;
t93 = -t67 * qJ(5) + t46;
t23 = -t142 * t93 + t33 * t85;
t30 = -t174 * t87 + t175 * t90;
t15 = -t45 * qJ(5) - t66 * qJD(5) - t30;
t31 = t174 * t90 + t175 * t87;
t92 = t44 * qJ(5) - t67 * qJD(5) + t31;
t6 = -t142 * t92 + t15 * t85;
t116 = t23 * t27 + t6 * t40;
t24 = t142 * t33 + t85 * t93;
t7 = t142 * t15 + t85 * t92;
t173 = t100 * t7 + t101 * t24 + t116;
t86 = sin(qJ(6));
t83 = t86 ^ 2;
t89 = cos(qJ(6));
t84 = t89 ^ 2;
t145 = t83 - t84;
t165 = t145 * qJD(6);
t81 = qJD(6) * t89;
t18 = t27 * t86 + t40 * t81;
t140 = qJD(6) * t86;
t170 = t40 * t140 - t27 * t89;
t76 = pkin(3) * t88 + qJ(2);
t51 = pkin(4) * t66 + t76;
t96 = pkin(5) * t100 + pkin(9) * t40 + t51;
t94 = t89 * t96;
t8 = -t86 * t24 + t94;
t9 = t89 * t24 + t86 * t96;
t117 = t8 * t86 - t89 * t9;
t138 = t91 * qJD(3);
t72 = pkin(3) * t138 + qJD(2);
t35 = t45 * pkin(4) + t72;
t95 = pkin(5) * t101 + pkin(9) * t27 + t35;
t2 = -qJD(6) * t94 + t24 * t140 - t89 * t7 - t86 * t95;
t3 = -t9 * qJD(6) - t86 * t7 + t89 * t95;
t162 = qJD(6) * t117 + t2 * t86 - t3 * t89;
t17 = t100 * t81 + t101 * t86;
t160 = t100 * t140 - t101 * t89;
t159 = 0.2e1 * qJD(2);
t157 = t23 * t6;
t156 = t23 * t81 + t6 * t86;
t155 = t23 * t54;
t151 = t40 * t86;
t150 = t40 * t89;
t149 = t66 * t45;
t148 = t67 * t44;
t147 = t86 * t89;
t146 = t52 * t81 + t54 * t86;
t144 = t83 + t84;
t137 = qJ(2) * qJD(3);
t135 = 0.2e1 * t166;
t132 = t75 * t140;
t131 = t75 * t81;
t130 = t86 * t81;
t129 = t88 * t138;
t128 = t158 * qJD(3);
t13 = t144 * t101;
t34 = t144 * t55;
t74 = pkin(4) * t85 + pkin(9);
t127 = t144 * t74;
t126 = 0.4e1 * t40 * t147;
t125 = t52 * t140 - t54 * t89;
t121 = t37 * t130;
t118 = -t8 * t89 - t86 * t9;
t115 = -t101 * t74 - t177;
t114 = -t100 ^ 2 - t37;
t53 = pkin(9) + t57;
t113 = t100 * t53 + t40 * t52;
t111 = t100 * t74 + t40 * t75;
t110 = t148 - t149;
t98 = -t101 * t53 + t112 - t178;
t97 = t30 * t66 - t31 * t67 + t44 * t46 - t45 * t47;
t1 = t118 * qJD(6) - t2 * t89 - t3 * t86;
t82 = qJ(2) * t159;
t71 = -0.2e1 * t130;
t70 = 0.2e1 * t130;
t65 = -0.2e1 * t165;
t20 = t23 * t140;
t12 = t147 * t27 - t40 * t165;
t10 = qJD(6) * t126 + t145 * t27;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t159, t82, -0.2e1 * t129, 0.2e1 * (t88 ^ 2 - t91 ^ 2) * qJD(3), 0, 0.2e1 * t129, 0, 0, 0.2e1 * qJD(2) * t88 + 0.2e1 * t91 * t137, 0.2e1 * qJD(2) * t91 - 0.2e1 * t137 * t88, 0, t82, -0.2e1 * t148, 0.2e1 * t44 * t66 - 0.2e1 * t45 * t67, 0, 0.2e1 * t149, 0, 0, 0.2e1 * t45 * t76 + 0.2e1 * t66 * t72, -0.2e1 * t44 * t76 + 0.2e1 * t67 * t72, 0.2e1 * t97, -0.2e1 * t30 * t47 + 0.2e1 * t31 * t46 + 0.2e1 * t72 * t76, 0.2e1 * t176, 0.2e1 * t100 * t27 + 0.2e1 * t101 * t40, 0, t135, 0, 0, 0.2e1 * t100 * t35 + 0.2e1 * t101 * t51, -0.2e1 * t27 * t51 - 0.2e1 * t35 * t40, -0.2e1 * t173, 0.2e1 * t24 * t7 + 0.2e1 * t35 * t51 + 0.2e1 * t157, 0.2e1 * t176 * t84 - 0.2e1 * t121, -t126 * t27 + 0.2e1 * t37 * t165, 0.2e1 * t100 * t170 - 0.2e1 * t101 * t150, 0.2e1 * t176 * t83 + 0.2e1 * t121, 0.2e1 * t100 * t18 + 0.2e1 * t101 * t151, t135, 0.2e1 * t100 * t3 + 0.2e1 * t101 * t8 - 0.2e1 * t6 * t151 - 0.2e1 * t18 * t23, 0.2e1 * t100 * t2 - 0.2e1 * t101 * t9 - 0.2e1 * t6 * t150 + 0.2e1 * t170 * t23, -0.2e1 * t118 * t27 - 0.2e1 * t162 * t40, -0.2e1 * t2 * t9 + 0.2e1 * t3 * t8 + 0.2e1 * t157; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t110, -t97, 0, 0, 0, 0, 0, 0, 0, 0, -t99, t173, 0, 0, 0, 0, 0, 0, t114 * t81 - t86 * t99, -t114 * t140 - t89 * t99, 0, t1 * t100 - t101 * t117 + t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t110, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t100 * t13 + 0.2e1 * t176; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t139, 0, -t138, 0, -t88 * t128, -t91 * t128, 0, 0, 0, 0, -t44, 0, -t45, 0, t31, t30, -t179 (-t30 * t87 + t31 * t90 + (-t46 * t87 + t47 * t90) * qJD(4)) * pkin(3), 0, 0, -t27, 0, -t101, 0, -t6, -t7, t181, t24 * t55 - t56 * t6 + t57 * t7 + t155, -t12, t10, t17, t12, -t160, 0, t20 + (-qJD(6) * t113 - t6) * t89 + t98 * t86, t113 * t140 + t89 * t98 + t156, t1, t1 * t53 - t117 * t55 + t6 * t52 + t155; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t139, -t138, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t45, 0, t179, 0, 0, 0, 0, 0, 0, -t27, -t101, 0, -t181, 0, 0, 0, 0, 0, 0, t170, t18, t13, t100 * t34 + t13 * t53 + t153 + t178; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t134, -0.2e1 * t133, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t54, -0.2e1 * t55, 0, -0.2e1 * t54 * t56 + 0.2e1 * t55 * t57, t70, t65, 0, t71, 0, 0, 0.2e1 * t125, 0.2e1 * t146, 0.2e1 * t34, 0.2e1 * t34 * t53 + 0.2e1 * t52 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, 0, -t45, 0, t31, t30, 0, 0, 0, 0, -t27, 0, -t101, 0, -t6, -t7, t182 (-t142 * t6 + t7 * t85) * pkin(4), -t12, t10, t17, t12, -t160, 0, t20 + t115 * t86 + (-qJD(6) * t111 - t6) * t89, t111 * t140 + t115 * t89 + t156, t1, t1 * t74 + t6 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t45, 0, 0, 0, 0, 0, 0, 0, 0, -t27, -t101, 0, -t182, 0, 0, 0, 0, 0, 0, t170, t18, t13, t101 * t127 + t177; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t134, -t133, 0, 0, 0, 0, 0, 0, 0, 0, -t54, -t55, 0 (-t142 * t54 + t55 * t85) * pkin(4), t70, t65, 0, t71, 0, 0, t125 + t132, t131 + t146, t34, t127 * t55 + t54 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, t65, 0, t71, 0, 0, 0.2e1 * t132, 0.2e1 * t131, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, -t27, 0, t35, 0, 0, 0, 0, 0, 0, -t160, -t17, t144 * t27, -t162; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t170, 0, t18, t101, t3, t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, t160, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, 0, -t140, 0, -t53 * t81 - t55 * t86, t140 * t53 - t55 * t89, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, 0, -t140, 0, -t74 * t81, t74 * t140, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t140, -t81, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t4;
