% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRPR7
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tau_reg [5x20]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPR7_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR7_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR7_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR7_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:19:53
% EndTime: 2019-12-31 18:19:57
% DurationCPUTime: 1.13s
% Computational Cost: add. (1608->226), mult. (3527->316), div. (0->0), fcn. (2478->14), ass. (0->136)
t94 = sin(pkin(8));
t75 = t94 * pkin(1) + pkin(6);
t151 = qJ(4) + t75;
t101 = cos(qJ(5));
t102 = cos(qJ(3));
t143 = qJD(1) * qJD(3);
t135 = t102 * t143;
t99 = sin(qJ(3));
t138 = t99 * t143;
t93 = sin(pkin(9));
t95 = cos(pkin(9));
t62 = t93 * t102 + t95 * t99;
t36 = t62 * qJDD(1) + t95 * t135 - t93 * t138;
t57 = t62 * qJD(1);
t98 = sin(qJ(5));
t44 = t98 * qJD(3) + t101 * t57;
t11 = t44 * qJD(5) - t101 * qJDD(3) + t98 * t36;
t149 = qJD(5) * t98;
t140 = t62 * t149;
t142 = t102 * qJDD(1);
t145 = t99 * qJDD(1);
t124 = t95 * t142 - t93 * t145;
t56 = t62 * qJD(3);
t34 = qJD(1) * t56 + qJDD(5) - t124;
t150 = qJD(1) * t99;
t152 = t95 * t102;
t55 = qJD(1) * t152 - t93 * t150;
t53 = qJD(5) - t55;
t120 = -t93 * t99 + t152;
t59 = t120 * qJD(3);
t173 = -(t34 * t62 + t53 * t59) * t101 + t53 * t140;
t96 = cos(pkin(8));
t170 = t96 * pkin(1);
t79 = t102 * pkin(3) + pkin(2);
t119 = -t79 - t170;
t109 = pkin(3) * t138 + t119 * qJDD(1) + qJDD(4);
t65 = t75 * qJDD(1);
t110 = qJ(4) * qJDD(1) + qJD(1) * qJD(4) + qJD(2) * qJD(3) + t65;
t128 = t151 * qJD(1);
t117 = t128 * qJD(3);
t86 = t102 * qJDD(2);
t16 = qJDD(3) * pkin(3) - t102 * t117 - t110 * t99 + t86;
t19 = (qJDD(2) - t117) * t99 + t110 * t102;
t3 = t95 * t16 - t93 * t19;
t1 = -qJDD(3) * pkin(4) - t3;
t46 = t99 * qJD(2) + t128 * t102;
t158 = t93 * t46;
t155 = qJD(3) * pkin(3);
t45 = t102 * qJD(2) - t128 * t99;
t40 = t45 + t155;
t17 = t95 * t40 - t158;
t13 = -qJD(3) * pkin(4) - t17;
t54 = t119 * qJD(1) + qJD(4);
t24 = -t55 * pkin(4) - t57 * pkin(7) + t54;
t4 = t93 * t16 + t95 * t19;
t133 = qJDD(3) * pkin(7) + qJD(5) * t24 + t4;
t130 = qJD(3) * t151;
t114 = -t99 * qJD(4) - t102 * t130;
t47 = t102 * qJD(4) - t99 * t130;
t23 = t93 * t114 + t95 * t47;
t28 = -pkin(4) * t120 - t62 * pkin(7) + t119;
t136 = t151 * t99;
t60 = t151 * t102;
t32 = -t93 * t136 + t95 * t60;
t172 = t1 * t62 + t13 * t59 - (qJD(5) * t28 + t23) * t53 + t133 * t120 - t32 * t34;
t90 = qJ(1) + pkin(8);
t81 = sin(t90);
t83 = cos(t90);
t126 = g(1) * t83 + g(2) * t81;
t74 = t93 * pkin(3) + pkin(7);
t89 = qJ(3) + pkin(9);
t80 = sin(t89);
t82 = cos(t89);
t171 = t126 * t80 - (pkin(3) * t150 + t57 * pkin(4) - t55 * pkin(7) + qJD(5) * t74) * t53 - g(3) * t82 - t1;
t146 = t101 * qJD(3);
t10 = qJD(5) * t146 + t98 * qJDD(3) + t101 * t36 - t57 * t149;
t169 = -t10 * t120 + t44 * t56;
t168 = g(3) * t102;
t167 = t10 * t98;
t166 = t13 * t62;
t165 = t28 * t34;
t42 = t98 * t57 - t146;
t164 = t42 * t53;
t163 = t44 * t53;
t162 = t44 * t57;
t161 = t57 * t42;
t160 = t81 * t98;
t159 = t83 * t98;
t38 = t95 * t46;
t157 = t98 * t34;
t18 = t93 * t40 + t38;
t91 = t99 ^ 2;
t156 = -t102 ^ 2 + t91;
t154 = t81 * t101;
t153 = t83 * t101;
t77 = -pkin(2) - t170;
t68 = qJD(1) * t77;
t148 = qJDD(2) - g(3);
t147 = qJD(5) * t101;
t141 = t99 * t155;
t14 = qJD(3) * pkin(7) + t18;
t35 = -qJD(3) * t57 + t124;
t8 = -t35 * pkin(4) - t36 * pkin(7) + t109;
t134 = qJD(5) * t14 - t8;
t129 = t101 * t53;
t125 = -g(1) * t81 + g(2) * t83;
t100 = sin(qJ(1));
t103 = cos(qJ(1));
t123 = g(1) * t100 - g(2) * t103;
t122 = t11 * t120 - t56 * t42;
t118 = t101 * t34 + (t55 * t98 - t149) * t53;
t116 = g(3) * t80 - t133;
t113 = -t68 * qJD(1) + t126 - t65;
t112 = 0.2e1 * t68 * qJD(3) - qJDD(3) * t75;
t21 = t95 * t45 - t158;
t111 = -t74 * t34 + (t13 + t21) * t53;
t104 = qJD(3) ^ 2;
t108 = -0.2e1 * qJDD(1) * t77 - t104 * t75 - t125;
t107 = (-t62 * t147 - t98 * t59) * t53 - t62 * t157;
t105 = qJD(1) ^ 2;
t97 = -qJ(4) - pkin(6);
t76 = -t95 * pkin(3) - pkin(4);
t64 = qJDD(3) * t102 - t104 * t99;
t63 = qJDD(3) * t99 + t104 * t102;
t51 = t82 * t153 + t160;
t50 = -t82 * t159 + t154;
t49 = -t82 * t154 + t159;
t48 = t82 * t160 + t153;
t31 = t95 * t136 + t93 * t60;
t26 = t56 * pkin(4) - t59 * pkin(7) + t141;
t22 = -t95 * t114 + t93 * t47;
t20 = t93 * t45 + t38;
t7 = t101 * t8;
t6 = t101 * t14 + t98 * t24;
t5 = t101 * t24 - t98 * t14;
t2 = [qJDD(1), t123, g(1) * t103 + g(2) * t100, (t123 + (t94 ^ 2 + t96 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t91 * qJDD(1) + 0.2e1 * t99 * t135, 0.2e1 * t99 * t142 - 0.2e1 * t156 * t143, t63, t64, 0, t108 * t102 + t112 * t99, t112 * t102 - t108 * t99, t120 * t4 - t17 * t59 - t18 * t56 + t22 * t57 + t23 * t55 - t3 * t62 + t31 * t36 + t32 * t35 - t126, t4 * t32 + t18 * t23 - t3 * t31 - t17 * t22 + t54 * t141 - g(1) * (-t100 * pkin(1) - t81 * t79 - t83 * t97) - g(2) * (t103 * pkin(1) + t83 * t79 - t81 * t97) + t109 * t119, -t44 * t140 + (t10 * t62 + t44 * t59) * t101, (-t101 * t42 - t44 * t98) * t59 + (-t167 - t101 * t11 + (-t101 * t44 + t42 * t98) * qJD(5)) * t62, t169 - t173, t107 + t122, -t120 * t34 + t53 * t56, -g(1) * t49 - g(2) * t51 + t31 * t11 + t22 * t42 + t5 * t56 - t7 * t120 + (t26 * t53 + t165 + (t120 * t14 - t32 * t53 + t166) * qJD(5)) * t101 + t172 * t98, -g(1) * t48 - g(2) * t50 + t31 * t10 + t22 * t44 - t6 * t56 + (-(-qJD(5) * t32 + t26) * t53 - t165 - t134 * t120 - qJD(5) * t166) * t98 + t172 * t101; 0, 0, 0, t148, 0, 0, 0, 0, 0, t64, -t63, -t120 * t36 + t62 * t35 + t59 * t55 + t56 * t57, t120 * t3 - t17 * t56 + t18 * t59 + t4 * t62 - g(3), 0, 0, 0, 0, 0, t107 - t122, t169 + t173; 0, 0, 0, 0, -t99 * t105 * t102, t156 * t105, t145, t142, qJDD(3), t113 * t99 - t168 + t86, t113 * t102 - t148 * t99, (t18 - t20) * t57 + (t17 - t21) * t55 + (t35 * t93 - t36 * t95) * pkin(3), t17 * t20 - t18 * t21 + (-t168 + t3 * t95 + t4 * t93 + (-qJD(1) * t54 + t126) * t99) * pkin(3), t129 * t44 + t167, (-t11 - t163) * t98 + (t10 - t164) * t101, t129 * t53 + t157 - t162, t118 + t161, -t53 * t57, t171 * t101 + t76 * t11 + t111 * t98 - t20 * t42 - t5 * t57, t76 * t10 + t111 * t101 - t171 * t98 - t20 * t44 + t6 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55 ^ 2 - t57 ^ 2, t17 * t57 - t18 * t55 + t109 + t125, 0, 0, 0, 0, 0, t118 - t161, -t101 * t53 ^ 2 - t157 - t162; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44 * t42, -t42 ^ 2 + t44 ^ 2, t10 + t164, -t11 + t163, t34, -g(1) * t50 + g(2) * t48 + t116 * t98 - t13 * t44 - t14 * t147 + t6 * t53 + t7, g(1) * t51 - g(2) * t49 + t101 * t116 + t13 * t42 + t134 * t98 + t5 * t53;];
tau_reg = t2;
