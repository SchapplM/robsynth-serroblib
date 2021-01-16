% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% tau_reg [5x18]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:56
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:55:47
% EndTime: 2021-01-15 16:55:52
% DurationCPUTime: 1.19s
% Computational Cost: add. (1100->217), mult. (2193->294), div. (0->0), fcn. (1410->10), ass. (0->133)
t72 = qJ(1) + pkin(7);
t66 = sin(t72);
t67 = cos(t72);
t142 = g(2) * t67 + g(3) * t66;
t78 = cos(pkin(7));
t62 = -t78 * pkin(1) - pkin(2);
t125 = qJDD(1) * t62;
t42 = qJDD(3) + t125;
t167 = -t142 - t42;
t75 = sin(pkin(8));
t119 = qJD(1) * qJD(4);
t82 = cos(qJ(4));
t109 = t82 * t119;
t80 = sin(qJ(4));
t89 = qJDD(1) * t80 + t109;
t166 = t89 * t75;
t77 = cos(pkin(8));
t128 = t77 * qJD(1);
t53 = -qJD(4) + t128;
t165 = qJD(4) + t53;
t120 = qJD(1) * qJD(3);
t76 = sin(pkin(7));
t57 = t76 * pkin(1) + qJ(3);
t39 = qJDD(1) * t57 + t120;
t164 = -t75 * (-qJ(5) - pkin(6)) + (t82 * pkin(4) + pkin(3)) * t77;
t147 = t77 * t80;
t26 = -t66 * t147 - t67 * t82;
t28 = t67 * t147 - t66 * t82;
t163 = -g(2) * t26 - g(3) * t28;
t158 = g(1) * t80;
t162 = t75 * t158 + t163;
t36 = -t77 * pkin(3) - t75 * pkin(6) + t62;
t25 = t36 * qJD(1) + qJD(3);
t44 = t57 * qJD(1);
t32 = t75 * qJD(2) + t77 * t44;
t108 = t82 * t25 - t80 * t32;
t137 = qJD(1) * t75;
t110 = qJ(5) * t137;
t6 = -t82 * t110 + t108;
t3 = -t53 * pkin(4) + t6;
t161 = -t6 + t3;
t159 = pkin(4) * t80;
t122 = t77 * qJDD(1);
t52 = -qJDD(4) + t122;
t155 = t52 * pkin(4);
t149 = t75 * t39;
t64 = t77 * qJDD(2);
t23 = -t64 + t149;
t154 = t23 * t75;
t153 = t52 * t77;
t151 = t66 * t80;
t70 = t75 ^ 2;
t84 = qJD(1) ^ 2;
t150 = t70 * t84;
t146 = t77 * t82;
t130 = qJD(4) * t82;
t133 = qJD(3) * t82;
t145 = t36 * t130 + t77 * t133;
t40 = t57 * t146;
t144 = t80 * t36 + t40;
t81 = sin(qJ(1));
t143 = t81 * pkin(1) + t66 * pkin(2);
t141 = t77 ^ 2 + t70;
t73 = t80 ^ 2;
t74 = t82 ^ 2;
t140 = -t73 - t74;
t139 = t73 - t74;
t138 = qJ(5) * t75;
t136 = qJD(1) * t80;
t135 = qJD(1) * t82;
t134 = qJD(3) * t77;
t132 = qJD(4) * t32;
t131 = qJD(4) * t80;
t129 = qJD(5) * t75;
t65 = t77 * qJD(2);
t15 = qJD(5) - t65 + (pkin(4) * t136 + t44) * t75;
t127 = qJD(5) + t15;
t123 = qJDD(1) * t82;
t121 = qJ(5) * qJDD(1);
t118 = qJD(1) * qJD(5);
t117 = t80 * t150;
t83 = cos(qJ(1));
t116 = t83 * pkin(1) + t67 * pkin(2) + t66 * qJ(3);
t115 = t82 * t138;
t114 = t75 * t136;
t113 = t57 * t131;
t112 = t53 * t131;
t31 = t75 * t44 - t65;
t111 = t31 * t137;
t22 = t36 * qJDD(1) + qJDD(3);
t24 = t75 * qJDD(2) + t77 * t39;
t107 = t25 * t130 - t32 * t131 + t80 * t22 + t82 * t24;
t35 = (t57 + t159) * t75;
t106 = qJD(1) * t35 + t15;
t105 = -qJD(4) * t25 - t24;
t104 = t52 - t122;
t103 = t52 + t122;
t102 = pkin(4) * t166 + qJDD(5) - t64;
t101 = qJD(4) * t114;
t100 = g(2) * t28 - g(3) * t26;
t27 = t66 * t146 - t67 * t80;
t29 = t67 * t146 + t151;
t99 = -g(2) * t29 - g(3) * t27;
t97 = -g(2) * t66 + g(3) * t67;
t96 = -g(2) * t83 - g(3) * t81;
t92 = -t80 * t25 - t82 * t32;
t7 = -t80 * t110 - t92;
t95 = t3 * t82 + t7 * t80;
t94 = t3 * t80 - t7 * t82;
t93 = t24 * t77 + t154;
t91 = t31 * t75 + t32 * t77;
t13 = t102 + t149;
t38 = (pkin(4) * t130 + qJD(3)) * t75;
t90 = qJD(1) * t38 + qJDD(1) * t35 + t13;
t18 = t82 * t22;
t88 = qJ(5) * t101 + t105 * t80 + t18;
t87 = -t53 ^ 2 - t150;
t86 = g(1) * t75 * t82 + g(2) * t27 - g(3) * t29 - t107;
t50 = t75 * t123;
t41 = t77 * t101;
t37 = t75 * t53 * t135;
t34 = t82 * t36;
t14 = -t80 * t138 + t144;
t12 = -t115 + t34 + (-t57 * t80 - pkin(4)) * t77;
t11 = t80 * t52 + t87 * t82;
t10 = -t82 * t52 + t87 * t80;
t9 = (t104 * t80 + (t53 - t128) * t130) * t75;
t8 = t41 + (t104 * t82 - t112) * t75;
t5 = -t80 * t134 - t82 * t129 + (-t40 + (-t36 + t138) * t80) * qJD(4);
t4 = -t80 * t129 + (-t57 * t147 - t115) * qJD(4) + t145;
t2 = (-t89 * qJ(5) - t80 * t118) * t75 + t107;
t1 = -t155 + (-t132 + (-t118 - t121) * t75) * t82 + t88;
t16 = [qJDD(1), t96, g(2) * t81 - g(3) * t83, (t96 + (t76 ^ 2 + t78 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), (-t125 + t167) * t77, t141 * t39 + t93 + t97, t42 * t62 - g(2) * t116 - g(3) * (-t67 * qJ(3) + t143) + t93 * t57 + t91 * qJD(3), (qJDD(1) * t74 - 0.2e1 * t80 * t109) * t70, 0.2e1 * (t139 * t119 - t80 * t123) * t70, t41 + (-t103 * t82 + t112) * t75, (t103 * t80 + (t53 + t128) * t130) * t75, t153, -t18 * t77 - t34 * t52 + ((qJD(1) * t70 + t53 * t77) * t57 + t91) * t130 + (-(-qJD(4) * t36 - t134) * t53 - t105 * t77 + t70 * t120 + t154 + (t70 * qJDD(1) + t153) * t57) * t80 + t99, (-t77 * t113 + t145) * t53 + t144 * t52 + t107 * t77 + (-t31 * t131 + t23 * t82) * t75 + (t57 * t123 + (-t113 + t133) * qJD(1)) * t70 + t100, -t1 * t77 - t12 * t52 - t5 * t53 + (t106 * t130 + t90 * t80) * t75 + t99, t14 * t52 + t2 * t77 + t4 * t53 + (-t106 * t131 + t90 * t82) * t75 + t100, ((-qJD(4) * t7 - qJDD(1) * t12 - t1 + (-qJD(4) * t14 - t5) * qJD(1)) * t82 + (qJD(4) * t3 - qJDD(1) * t14 - t2 + (qJD(4) * t12 - t4) * qJD(1)) * t80 - t142) * t75, t2 * t14 + t7 * t4 + t1 * t12 + t3 * t5 + t13 * t35 + t15 * t38 - g(2) * (pkin(4) * t151 + t116) - g(3) * (t164 * t66 + t143) + (-g(2) * t164 - g(3) * (-qJ(3) - t159)) * t67; 0, 0, 0, qJDD(2) - g(1), 0, 0, -t23 * t77 + t24 * t75 - g(1), 0, 0, 0, 0, 0, t9, t8, t9, t8, 0, -t13 * t77 - g(1) + (-t95 * qJD(4) - t1 * t80 + t2 * t82) * t75; 0, 0, 0, 0, -t122, -t141 * t84, -t91 * qJD(1) - t167, 0, 0, 0, 0, 0, t10, t11, t10, t11, t140 * t75 * qJDD(1), t1 * t82 + t2 * t80 - t94 * qJD(4) + (-t15 * t75 + t94 * t77) * qJD(1) + t142; 0, 0, 0, 0, 0, 0, 0, t82 * t117, -t139 * t150, -t114 * t165 + t50, -t37 - t166, -t52, -t82 * t111 + t165 * t92 - t80 * t24 + t162 + t18, -t108 * t53 + t80 * t111 + t86, -0.2e1 * t155 - t7 * t53 + (-pkin(4) * t117 - t132 + (-t127 * qJD(1) - t121) * t75) * t82 + t88 + t162, -t74 * pkin(4) * t150 - t6 * t53 + (t80 * t121 + (qJ(5) * t130 + t127 * t80) * qJD(1)) * t75 + t86, (-pkin(4) * t123 + (pkin(4) * qJD(4) - t161) * t136) * t75, t161 * t7 + (t1 + (-t15 * t135 + t158) * t75 + t163) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37 + t166, t50 + (-qJD(4) + t53) * t114, t140 * t150, g(1) * t77 + (t95 * qJD(1) + t39 + t97) * t75 + t102;];
tau_reg = t16;
