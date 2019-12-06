% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRPRP4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP4_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:36:13
% EndTime: 2019-12-05 15:36:16
% DurationCPUTime: 1.02s
% Computational Cost: add. (1099->206), mult. (2229->254), div. (0->0), fcn. (1591->10), ass. (0->134)
t128 = qJD(1) * qJD(2);
t89 = sin(qJ(2));
t91 = cos(qJ(2));
t172 = t89 * qJDD(1) + t91 * t128;
t143 = qJDD(2) * pkin(2);
t78 = t91 * qJDD(1);
t40 = -t89 * t128 + t143 + t78;
t84 = sin(pkin(8));
t86 = cos(pkin(8));
t15 = t172 * t86 + t84 * t40;
t174 = qJDD(2) * pkin(6) + qJD(3) * qJD(4) + t15;
t148 = qJD(1) * t89;
t138 = t91 * qJD(1);
t59 = qJD(2) * pkin(2) + t138;
t27 = t86 * t148 + t84 * t59;
t22 = qJD(2) * pkin(6) + t27;
t88 = sin(qJ(4));
t154 = t88 * t22;
t90 = cos(qJ(4));
t18 = t90 * qJD(3) - t154;
t173 = qJD(5) - t18;
t110 = t90 * pkin(4) + t88 * qJ(5);
t105 = pkin(3) + t110;
t160 = t86 * pkin(2);
t39 = -t105 - t160;
t14 = -t172 * t84 + t86 * t40;
t109 = pkin(4) * t88 - qJ(5) * t90;
t34 = t109 * qJD(4) - t88 * qJD(5);
t147 = qJD(2) * t34;
t8 = -t105 * qJDD(2) - t14 + t147;
t171 = -qJDD(2) * t39 - t8;
t85 = sin(pkin(7));
t87 = cos(pkin(7));
t112 = g(1) * t87 + g(2) * t85;
t81 = qJ(2) + pkin(8);
t73 = sin(t81);
t170 = t112 * t73;
t13 = -qJD(4) * pkin(4) + t173;
t152 = t90 * t22;
t19 = t88 * qJD(3) + t152;
t16 = qJD(4) * qJ(5) + t19;
t43 = t84 * t89 - t86 * t91;
t146 = qJD(2) * t43;
t44 = t84 * t91 + t86 * t89;
t169 = -qJD(2) * t146 + t44 * qJDD(2);
t142 = qJDD(4) * pkin(4);
t168 = qJDD(5) - t142;
t162 = g(3) * t73;
t74 = cos(t81);
t167 = t112 * t74 + t162;
t166 = pkin(2) * t89;
t165 = pkin(6) * t74;
t161 = g(3) * t74;
t69 = t84 * pkin(2) + pkin(6);
t92 = qJD(4) ^ 2;
t159 = t69 * t92;
t158 = t85 * t88;
t157 = t85 * t90;
t156 = t87 * t88;
t155 = t87 * t90;
t153 = t88 * t90;
t36 = t44 * qJD(1);
t151 = t34 - t36;
t82 = t88 ^ 2;
t83 = t90 ^ 2;
t150 = -t82 + t83;
t149 = t82 + t83;
t145 = qJD(2) * t88;
t144 = qJD(4) * t88;
t141 = t36 * qJD(2);
t62 = t84 * t148;
t38 = t86 * t138 - t62;
t139 = t38 * qJD(2);
t137 = qJDD(1) - g(3);
t134 = qJDD(4) * t69;
t132 = t82 * qJDD(2);
t131 = t83 * qJDD(2);
t129 = t90 * qJDD(2);
t127 = qJD(2) * qJD(4);
t125 = qJDD(4) * qJ(5);
t124 = t88 * qJDD(3) + t174 * t90;
t123 = t91 * pkin(2) + t74 * pkin(3) + t73 * pkin(6);
t122 = -g(1) * t85 + g(2) * t87;
t120 = t18 + t154;
t26 = t86 * t59 - t62;
t119 = qJD(4) * t152 - t90 * qJDD(3) + t174 * t88;
t118 = t90 * t141 + t38 * t144 + (g(1) * t155 + g(2) * t157) * t73;
t17 = -t105 * qJD(2) - t26;
t117 = qJD(2) * t39 + t17;
t21 = -qJD(2) * pkin(3) - t26;
t70 = -pkin(3) - t160;
t116 = qJD(2) * t70 + t21;
t114 = t127 * t153;
t113 = -pkin(3) * t73 - t166;
t111 = t159 + t161;
t108 = t13 * t88 + t16 * t90;
t107 = t18 * t88 - t19 * t90;
t35 = t44 * qJD(2);
t106 = -t35 * qJD(2) - t43 * qJDD(2);
t31 = t74 * t157 - t156;
t33 = t74 * t155 + t158;
t104 = g(1) * t33 + g(2) * t31 - t124;
t103 = t44 * t92 - t106;
t11 = -qJDD(2) * pkin(3) - t14;
t102 = qJDD(2) * t70 + t11 + t111;
t101 = 0.2e1 * t146 * qJD(4) - qJDD(4) * t44;
t100 = -t161 + t170;
t99 = -g(3) * t91 + t112 * t89;
t30 = t74 * t158 + t155;
t32 = t74 * t156 - t157;
t98 = g(1) * t32 + g(2) * t30 + t88 * t162 - t119;
t97 = t19 * qJD(4) + t98;
t3 = t125 + (qJD(5) - t154) * qJD(4) + t124;
t4 = t119 + t168;
t96 = t3 * t90 + t4 * t88 + (t13 * t90 - t16 * t88) * qJD(4);
t5 = -t22 * t144 + t124;
t95 = t5 * t90 + t119 * t88 + (-t18 * t90 - t19 * t88) * qJD(4);
t94 = -t149 * t139 - t167 + (t131 + t132) * t69;
t93 = qJD(2) ^ 2;
t76 = t88 * qJDD(2);
t63 = t93 * t153;
t55 = t87 * t165;
t54 = t85 * t165;
t53 = t150 * t93;
t52 = qJDD(4) * t90 - t92 * t88;
t51 = qJDD(4) * t88 + t92 * t90;
t45 = t109 * qJD(2);
t42 = -0.2e1 * t114 + t131;
t41 = 0.2e1 * t114 + t132;
t25 = t150 * t127 + t88 * t129;
t7 = t169 * t149;
t2 = t101 * t88 - t103 * t90;
t1 = t101 * t90 + t103 * t88;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t137, 0, 0, 0, 0, 0, 0, t91 * qJDD(2) - t93 * t89, -qJDD(2) * t89 - t93 * t91, 0, -g(3) + (t89 ^ 2 + t91 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, t106, -t169, 0, -t14 * t43 - t146 * t27 + t15 * t44 - t26 * t35 - g(3), 0, 0, 0, 0, 0, 0, t2, t1, t7, t107 * t146 + t11 * t43 + t21 * t35 + t95 * t44 - g(3), 0, 0, 0, 0, 0, 0, t2, t7, -t1, -t108 * t146 + t17 * t35 + t8 * t43 + t96 * t44 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t78 + t99, t112 * t91 - t137 * t89, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t86 * t143 + t100 + t14 + t141, -t84 * t143 + t139 - t15 + t167, 0, t26 * t36 - t27 * t38 + (t14 * t86 + t15 * t84 + t99) * pkin(2), t41, 0.2e1 * t25, t51, t42, t52, 0, (t116 * qJD(4) - t134) * t88 - t102 * t90 + t118, (-t134 + (t116 + t38) * qJD(4)) * t90 + (t102 - t141 - t170) * t88, t94 + t95, t11 * t70 - t21 * t36 - g(1) * (t113 * t87 + t55) - g(2) * (t113 * t85 + t54) - g(3) * t123 + t107 * t38 + t95 * t69, t41, t51, -0.2e1 * t25, 0, -t52, t42, (t117 * qJD(4) - t134) * t88 + (-t111 - t147 + t171) * t90 + t118, t94 + t96, (t134 + (-t117 - t38) * qJD(4)) * t90 + (-t151 * qJD(2) + t100 - t159 + t171) * t88, t8 * t39 - g(1) * (-t87 * t166 + t55) - g(2) * (-t85 * t166 + t54) - g(3) * (t110 * t74 + t123) - t108 * t38 + t151 * t17 + t96 * t69 + t105 * t170; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3) + t122, 0, 0, 0, 0, 0, 0, t52, -t51, 0, -t107 * qJD(4) - t119 * t90 + t5 * t88 + t122, 0, 0, 0, 0, 0, 0, t52, 0, t51, t108 * qJD(4) + t3 * t88 - t4 * t90 + t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, -t53, t76, t63, t129, qJDD(4), -t21 * t145 + t97, (-qJD(2) * t21 + t162) * t90 + t120 * qJD(4) + t104, 0, 0, -t63, t76, t53, qJDD(4), -t129, t63, 0.2e1 * t142 - qJDD(5) + (-t17 * t88 + t45 * t90) * qJD(2) + t97, -t109 * qJDD(2), -t90 * t162 + 0.2e1 * t125 + (t17 * t90 + t45 * t88) * qJD(2) + (0.2e1 * qJD(5) - t120) * qJD(4) - t104, t3 * qJ(5) - t4 * pkin(4) - t17 * t45 - t13 * t19 - g(1) * (-t32 * pkin(4) + t33 * qJ(5)) - g(2) * (-t30 * pkin(4) + t31 * qJ(5)) + t109 * t162 + t173 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(4) - t63, t76, -t82 * t93 - t92, -t16 * qJD(4) + t17 * t145 + t168 - t98;];
tau_reg = t6;
