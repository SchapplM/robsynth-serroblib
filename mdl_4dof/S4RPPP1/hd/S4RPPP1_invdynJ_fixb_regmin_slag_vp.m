% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
% 
% Output:
% tau_reg [4x15]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:46
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau_reg = S4RPPP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:45:37
% EndTime: 2018-11-14 13:45:38
% DurationCPUTime: 0.50s
% Computational Cost: add. (533->175), mult. (1283->226), div. (0->0), fcn. (939->10), ass. (0->131)
t78 = cos(pkin(4));
t150 = pkin(1) * t78;
t79 = sin(qJ(1));
t149 = g(1) * t79;
t75 = sin(pkin(6));
t71 = t75 ^ 2;
t76 = sin(pkin(4));
t72 = t76 ^ 2;
t148 = t71 * t72;
t81 = qJD(1) ^ 2;
t147 = t72 * t81;
t146 = t75 * t76;
t77 = cos(pkin(6));
t145 = t76 * t77;
t80 = cos(qJ(1));
t144 = t76 * t80;
t143 = t78 * t81;
t130 = qJDD(1) * t76;
t113 = qJ(2) * t130;
t128 = qJD(1) * qJD(2);
t111 = t76 * t128;
t46 = t75 * t111;
t142 = t75 * t113 + t46;
t135 = qJD(1) * t76;
t117 = qJ(2) * t135;
t134 = qJD(1) * t78;
t119 = t75 * t134;
t31 = pkin(1) * t119 + t77 * t117;
t140 = qJ(2) * t76;
t36 = t77 * t140 + t75 * t150;
t141 = t80 * pkin(1) + t79 * t140;
t139 = t78 * qJ(3);
t94 = pkin(3) * t145 + t139;
t8 = t94 * qJD(1) + qJD(4) + t31;
t138 = -qJD(4) - t8;
t137 = pkin(1) * qJDD(1);
t133 = qJD(2) * t76;
t41 = t78 * qJD(3) + t77 * t133;
t136 = qJD(1) * t41;
t49 = t75 * t117;
t132 = qJD(3) + t49;
t131 = qJDD(1) * t75;
t129 = t78 * qJDD(1);
t127 = qJD(1) * qJD(3);
t70 = pkin(4) - pkin(6);
t63 = cos(t70) / 0.2e1;
t69 = pkin(4) + pkin(6);
t67 = cos(t69);
t126 = t63 + t67 / 0.2e1;
t62 = -sin(t70) / 0.2e1;
t66 = sin(t69);
t125 = t66 / 0.2e1 + t62;
t124 = pkin(3) * t146;
t123 = t77 * t150;
t122 = t77 * t143;
t121 = pkin(1) * t129;
t48 = t77 * t111;
t12 = t77 * t113 + t75 * t121 + t48;
t120 = t75 * t135;
t118 = -pkin(1) * t77 - pkin(2);
t116 = t77 * t130;
t115 = qJDD(3) + t142;
t114 = -t79 * pkin(1) + t80 * t140;
t112 = t72 * t128;
t110 = t75 * t127;
t109 = -qJ(3) * t75 - pkin(1);
t108 = t75 * t77 * t147;
t104 = t118 * t78;
t56 = t75 * t140;
t25 = t56 + t104;
t98 = t118 * qJDD(1);
t6 = t78 * t98 + t115;
t107 = qJDD(1) * t25 + t6;
t95 = -pkin(2) * t77 + t109;
t26 = t95 * t76;
t86 = t95 * qJDD(1);
t7 = qJDD(2) + (t86 - t110) * t76;
t106 = qJDD(1) * t26 + t7;
t105 = -qJ(4) + t118;
t4 = -qJ(3) * t129 - t78 * t127 - t12;
t20 = -t80 * t126 + t79 * t75;
t22 = t79 * t126 + t80 * t75;
t103 = g(1) * t20 - g(2) * t22;
t21 = t80 * t125 + t79 * t77;
t23 = -t79 * t125 + t80 * t77;
t102 = g(1) * t21 - g(2) * t23;
t101 = g(1) * t80 + g(2) * t79;
t100 = g(2) * t144 - g(3) * t78 + qJDD(2);
t99 = (qJD(1) * t123 - t49) * t75 - t31 * t77;
t97 = t23 * pkin(2) + t22 * qJ(3) + t141;
t96 = -qJD(3) * t75 - qJD(4) * t77;
t93 = t105 * qJDD(1);
t55 = -pkin(1) * t130 + qJDD(2);
t92 = t72 * t137 - t55 * t76;
t51 = qJDD(1) * t124;
t1 = t51 + (-qJD(1) * qJD(4) + t93) * t78 + t115;
t83 = t105 * t78 + t124;
t10 = t56 + t83;
t40 = -t78 * qJD(4) + t75 * t133;
t91 = qJD(1) * t40 + qJDD(1) * t10 + t1;
t13 = t94 + t36;
t3 = pkin(3) * t116 + qJDD(4) - t4;
t90 = qJDD(1) * t13 + t136 + t3;
t24 = -t139 - t36;
t89 = -qJDD(1) * t24 + t136 - t4;
t88 = -t21 * pkin(2) - t20 * qJ(3) + t114;
t87 = (-pkin(2) - qJ(4)) * t77 + t109;
t14 = t87 * t76;
t84 = t87 * qJDD(1);
t2 = qJDD(2) + (t96 * qJD(1) + t84) * t76;
t34 = t96 * t76;
t85 = (-qJD(1) * t34 - qJDD(1) * t14 - t2) * t76;
t82 = -g(1) * t22 - g(2) * t20 - g(3) * (-t66 / 0.2e1 + t62) + t115;
t74 = t78 ^ 2;
t73 = t77 ^ 2;
t44 = t71 * t112;
t43 = t143 * t146;
t39 = (-t74 - t148) * t81;
t35 = -t56 + t123;
t33 = (-t71 - t73) * t147;
t32 = t108 + t129;
t29 = -t43 + t116;
t28 = (-t122 + t131) * t76;
t27 = (t122 + t131) * t76;
t17 = -qJ(3) * t134 - t31;
t16 = qJD(1) * t26 + qJD(2);
t15 = qJD(1) * t104 + t132;
t11 = t77 * t121 - t142;
t9 = qJD(1) * t14 + qJD(2);
t5 = t83 * qJD(1) + t132;
t18 = [qJDD(1), -g(2) * t80 + t149, t101, t92 * t77 + (qJDD(1) * t35 + t11 - t46) * t78 + t102, -t92 * t75 + (-qJDD(1) * t36 - t12 - t48) * t78 - t103, t73 * t112 + t44 + (-t11 * t75 + t12 * t77 + (-t35 * t75 + t36 * t77) * qJDD(1) - t101) * t76, t12 * t36 + t11 * t35 - g(1) * t114 - g(2) * t141 + (-t55 * pkin(1) - t99 * qJD(2)) * t76, t44 + (t107 * t75 + t89 * t77 - t101) * t76, -t72 * t77 * t110 + t107 * t78 + (qJD(2) * t119 + t106 * t77) * t76 - t102, -t106 * t146 + t127 * t148 + t89 * t78 + t103, t7 * t26 + t4 * t24 - t17 * t41 + t6 * t25 - g(1) * t88 - g(2) * t97 + (qJD(2) * t15 - qJD(3) * t16) * t146 (t91 * t75 + t90 * t77 - t101) * t76, t75 * t85 + t90 * t78 + t103, t77 * t85 - t91 * t78 + t102, t2 * t14 + t9 * t34 + t1 * t10 + t5 * t40 + t3 * t13 + t8 * t41 - g(1) * (pkin(3) * t144 - t21 * qJ(4) + t88) - g(2) * (t79 * t76 * pkin(3) + t23 * qJ(4) + t97); 0, 0, 0, -t29, t27, t33 (t99 * qJD(1) - t137 - t149) * t76 + t100, t33, t29, -t27 (-t149 + t86 + (t17 * t77 + (-qJD(3) - t15) * t75) * qJD(1)) * t76 + t100, t33, -t27, -t29 (-t149 + t84 + (t138 * t77 + (-qJD(3) - t5) * t75) * qJD(1)) * t76 + t100; 0, 0, 0, 0, 0, 0, 0, t28, t32, t39, t16 * t120 + (qJD(1) * t17 + t98) * t78 + t82, t28, t39, -t32, t9 * t120 + t51 + (t138 * qJD(1) + t93) * t78 + t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43 + t116, -t108 + t129 (-t72 * t73 - t74) * t81, -g(1) * t23 - g(2) * t21 - g(3) * (t63 - t67 / 0.2e1) + (t9 * t145 + t5 * t78) * qJD(1) + t3;];
tau_reg  = t18;
