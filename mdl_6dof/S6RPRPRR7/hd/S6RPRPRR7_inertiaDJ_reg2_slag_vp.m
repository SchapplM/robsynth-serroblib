% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RPRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRPRR7_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR7_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_inertiaDJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:56:04
% EndTime: 2019-03-09 03:56:12
% DurationCPUTime: 2.70s
% Computational Cost: add. (4932->191), mult. (9377->319), div. (0->0), fcn. (9484->8), ass. (0->121)
t141 = sin(pkin(10));
t142 = cos(pkin(10));
t87 = sin(qJ(3));
t89 = cos(qJ(3));
t101 = t141 * t87 - t142 * t89;
t146 = t141 * t89 + t142 * t87;
t158 = sin(qJ(5));
t159 = cos(qJ(5));
t194 = t159 * t101 + t158 * t146;
t115 = t142 * pkin(3) + pkin(4);
t126 = t141 * pkin(3);
t60 = t158 * t115 + t159 * t126;
t49 = t60 * qJD(5);
t153 = t194 * t49;
t167 = t159 * t146;
t175 = -t158 * t101 + t167;
t59 = t159 * t115 - t158 * t126;
t48 = t59 * qJD(5);
t112 = -t175 * t48 - t153;
t124 = qJD(5) * t158;
t147 = t146 * qJD(3);
t166 = t101 * qJD(3);
t171 = qJD(5) * t167 - t101 * t124 + t159 * t147 - t158 * t166;
t23 = qJD(5) * t194 + t158 * t147 + t159 * t166;
t198 = t171 * t59 + t23 * t60 + t112;
t90 = -pkin(1) - pkin(7);
t143 = qJ(4) - t90;
t68 = t143 * t87;
t69 = t143 * t89;
t42 = -t141 * t69 - t142 * t68;
t34 = -t146 * pkin(8) + t42;
t41 = t141 * t68 - t142 * t69;
t105 = pkin(8) * t101 + t41;
t96 = t159 * t105;
t19 = t158 * t34 - t96;
t20 = t158 * t105 + t159 * t34;
t138 = t87 * qJD(3);
t102 = -t89 * qJD(4) + t143 * t138;
t103 = -qJD(3) * t69 - t87 * qJD(4);
t33 = t141 * t102 + t142 * t103;
t28 = pkin(8) * t166 + t33;
t32 = t142 * t102 - t141 * t103;
t92 = -t147 * pkin(8) - t32;
t6 = t20 * qJD(5) + t158 * t28 + t159 * t92;
t114 = t171 * t19 + t194 * t6;
t91 = -qJD(5) * t96 + t34 * t124 + t158 * t92 - t159 * t28;
t197 = -t175 * t91 - t20 * t23 + t114;
t195 = t194 * t171;
t196 = t175 * t23;
t193 = t195 - t196;
t36 = t194 ^ 2;
t86 = sin(qJ(6));
t139 = qJD(6) * t86;
t88 = cos(qJ(6));
t106 = t139 * t175 + t23 * t88;
t82 = qJD(6) * t88;
t191 = -t175 * t82 + t23 * t86;
t108 = -t171 * t86 - t194 * t82;
t156 = t171 * t88;
t107 = -t139 * t194 + t156;
t189 = pkin(5) * t171;
t188 = (t141 * t166 + t142 * t147) * pkin(3);
t55 = -pkin(5) - t59;
t186 = t171 * t55;
t128 = t101 * t147;
t176 = t146 * t166;
t177 = -0.2e1 * t128 + 0.2e1 * t176;
t84 = t86 ^ 2;
t85 = t88 ^ 2;
t145 = t84 - t85;
t165 = t145 * qJD(6);
t174 = t101 * t32 - t33 * t146 + t41 * t147 + t166 * t42;
t81 = t87 * pkin(3) + qJ(2);
t47 = t146 * pkin(4) + t81;
t95 = pkin(5) * t175 + pkin(9) * t194 + t47;
t170 = -qJD(6) * t95 + t91;
t8 = -t86 * t20 + t88 * t95;
t9 = t88 * t20 + t86 * t95;
t120 = t8 * t86 - t88 * t9;
t137 = t89 * qJD(3);
t74 = pkin(3) * t137 + qJD(2);
t43 = -pkin(4) * t166 + t74;
t93 = -pkin(5) * t23 + pkin(9) * t171 + t43;
t2 = t20 * t139 + t170 * t88 - t86 * t93;
t3 = t170 * t86 - t20 * t82 + t88 * t93;
t164 = qJD(6) * t120 + t2 * t86 - t3 * t88;
t162 = 0.2e1 * qJD(2);
t161 = t19 * t6;
t160 = t19 * t82 + t6 * t86;
t157 = t19 * t49;
t151 = t194 * t86;
t150 = t194 * t88;
t149 = t49 * t86 + t55 * t82;
t144 = t84 + t85;
t136 = qJ(2) * qJD(3);
t135 = -0.2e1 * t196;
t134 = t86 * t156;
t133 = pkin(5) * t139;
t132 = pkin(5) * t82;
t131 = t86 * t82;
t130 = t87 * t137;
t12 = t144 * t23;
t31 = t144 * t48;
t56 = pkin(9) + t60;
t127 = t144 * t56;
t125 = t55 * t139 - t49 * t88;
t123 = t36 * t131;
t121 = t8 * t88 + t86 * t9;
t113 = -t175 ^ 2 - t36;
t111 = t175 * t56 + t194 * t55;
t100 = -0.2e1 * t193;
t99 = t23 * t56 + t112 - t186;
t1 = -t121 * qJD(6) - t2 * t88 - t3 * t86;
t83 = qJ(2) * t162;
t73 = -0.2e1 * t131;
t72 = 0.2e1 * t131;
t67 = -0.2e1 * t165;
t17 = t19 * t139;
t11 = -t165 * t194 + t134;
t7 = 0.4e1 * t131 * t194 + t145 * t171;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t162, t83, -0.2e1 * t130, 0.2e1 * (t87 ^ 2 - t89 ^ 2) * qJD(3), 0, 0.2e1 * t130, 0, 0, 0.2e1 * qJD(2) * t87 + 0.2e1 * t89 * t136, 0.2e1 * qJD(2) * t89 - 0.2e1 * t87 * t136, 0, t83, 0.2e1 * t128, -0.2e1 * t101 * t166 + 0.2e1 * t147 * t146, 0, -0.2e1 * t176, 0, 0, 0.2e1 * t74 * t146 - 0.2e1 * t166 * t81, -0.2e1 * t101 * t74 - 0.2e1 * t81 * t147, 0.2e1 * t174, 0.2e1 * t32 * t41 + 0.2e1 * t33 * t42 + 0.2e1 * t74 * t81, 0.2e1 * t195, 0.2e1 * t171 * t175 - 0.2e1 * t194 * t23, 0, t135, 0, 0, 0.2e1 * t175 * t43 - 0.2e1 * t23 * t47, -0.2e1 * t171 * t47 - 0.2e1 * t194 * t43, -0.2e1 * t197, -0.2e1 * t20 * t91 + 0.2e1 * t47 * t43 + 0.2e1 * t161, 0.2e1 * t195 * t85 - 0.2e1 * t123, -0.4e1 * t134 * t194 + 0.2e1 * t165 * t36, -0.2e1 * t107 * t175 + 0.2e1 * t150 * t23, 0.2e1 * t195 * t84 + 0.2e1 * t123, -0.2e1 * t108 * t175 - 0.2e1 * t151 * t23, t135, 0.2e1 * t108 * t19 - 0.2e1 * t6 * t151 + 0.2e1 * t175 * t3 - 0.2e1 * t23 * t8, -0.2e1 * t107 * t19 - 0.2e1 * t6 * t150 + 0.2e1 * t175 * t2 + 0.2e1 * t23 * t9, 0.2e1 * t121 * t171 - 0.2e1 * t164 * t194, -0.2e1 * t2 * t9 + 0.2e1 * t3 * t8 + 0.2e1 * t161; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t177, -t174, 0, 0, 0, 0, 0, 0, 0, 0, t100, t197, 0, 0, 0, 0, 0, 0, t100 * t86 + t113 * t82, t100 * t88 - t113 * t139, 0, t1 * t175 + t120 * t23 + t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t177, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t193, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t12 * t175 + 0.2e1 * t195; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t138, 0, -t137, 0, -t90 * t138, -t90 * t137, 0, 0, 0, 0, -t147, 0, t166, 0, t32, -t33, t188 (t141 * t33 + t142 * t32) * pkin(3), 0, 0, -t171, 0, t23, 0, -t6, t91, t198, t20 * t48 - t6 * t59 - t60 * t91 + t157, -t11, t7, -t191, t11, -t106, 0, t17 + (-t111 * qJD(6) - t6) * t88 + t99 * t86, t111 * t139 + t88 * t99 + t160, t1, t1 * t56 - t120 * t48 + t55 * t6 + t157; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t138, -t137, 0, 0, 0, 0, 0, 0, 0, 0, -t147, t166, 0, -t188, 0, 0, 0, 0, 0, 0, -t171, t23, 0, -t198, 0, 0, 0, 0, 0, 0, -t107, -t108, -t12, -t127 * t23 + t175 * t31 + t153 + t186; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t49, -0.2e1 * t48, 0, 0.2e1 * t48 * t60 - 0.2e1 * t49 * t59, t72, t67, 0, t73, 0, 0, 0.2e1 * t125, 0.2e1 * t149, 0.2e1 * t31, 0.2e1 * t127 * t48 + 0.2e1 * t55 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t166, -t147, 0, t74, 0, 0, 0, 0, 0, 0, -t23, -t171, 0, t43, 0, 0, 0, 0, 0, 0, -t106, t191, t144 * t171, -t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t171, 0, t23, 0, -t6, t91, 0, 0, -t11, t7, -t191, t11, -t106, 0, t17 + (pkin(9) * t23 + t189) * t86 + (-t6 + (pkin(5) * t194 - pkin(9) * t175) * qJD(6)) * t88, pkin(5) * t107 + pkin(9) * t106 + t160, t1, -t6 * pkin(5) + pkin(9) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t171, t23, 0, 0, 0, 0, 0, 0, 0, 0, -t107, -t108, -t12, -pkin(9) * t12 - t189; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t48, 0, 0, t72, t67, 0, t73, 0, 0, t125 - t133, -t132 + t149, t31, -t49 * pkin(5) + pkin(9) * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, t67, 0, t73, 0, 0, -0.2e1 * t133, -0.2e1 * t132, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t107, 0, -t108, -t23, t3, t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t191, t106, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, 0, -t139, 0, -t48 * t86 - t56 * t82, t56 * t139 - t48 * t88, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t139, -t82, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, 0, -t139, 0, -pkin(9) * t82, pkin(9) * t139, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t4;