% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPPRR6
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% tau_reg [5x22]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRR6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR6_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR6_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR6_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:58:09
% EndTime: 2019-12-31 17:58:12
% DurationCPUTime: 0.99s
% Computational Cost: add. (1356->209), mult. (2941->278), div. (0->0), fcn. (2231->14), ass. (0->123)
t87 = sin(pkin(9));
t92 = sin(qJ(4));
t135 = t92 * t87;
t123 = qJD(1) * t135;
t89 = cos(pkin(9));
t133 = qJD(1) * t89;
t95 = cos(qJ(4));
t67 = t95 * t133;
t50 = t67 - t123;
t47 = qJD(5) - t50;
t86 = qJ(1) + pkin(8);
t77 = sin(t86);
t79 = cos(t86);
t116 = g(1) * t79 + g(2) * t77;
t85 = pkin(9) + qJ(4);
t76 = sin(t85);
t78 = cos(t85);
t101 = -g(3) * t78 + t116 * t76;
t88 = sin(pkin(8));
t68 = t88 * pkin(1) + qJ(3);
t56 = qJD(1) * qJD(3) + t68 * qJDD(1);
t73 = t89 * qJDD(2);
t31 = t73 + (-pkin(6) * qJDD(1) - t56) * t87;
t127 = t89 * qJDD(1);
t39 = t87 * qJDD(2) + t89 * t56;
t32 = pkin(6) * t127 + t39;
t111 = -t95 * t31 + t92 * t32;
t62 = t68 * qJD(1);
t75 = t89 * qJD(2);
t33 = t75 + (-pkin(6) * qJD(1) - t62) * t87;
t45 = t87 * qJD(2) + t89 * t62;
t34 = pkin(6) * t133 + t45;
t13 = t92 * t33 + t95 * t34;
t2 = -qJDD(4) * pkin(4) + t13 * qJD(4) + t111;
t59 = t95 * t87 + t92 * t89;
t51 = t59 * qJD(1);
t156 = t101 - (t51 * pkin(4) + t47 * pkin(7)) * t47 - t2;
t115 = g(1) * t77 - g(2) * t79;
t90 = cos(pkin(8));
t70 = -t90 * pkin(1) - pkin(2);
t129 = qJDD(1) * t70;
t60 = qJDD(3) + t129;
t155 = t115 - t60;
t91 = sin(qJ(5));
t132 = qJD(5) * t91;
t58 = -t95 * t89 + t135;
t52 = t58 * qJD(4);
t94 = cos(qJ(5));
t104 = t59 * t132 + t94 * t52;
t128 = t87 * qJDD(1);
t113 = -t95 * t127 + t92 * t128;
t53 = t59 * qJD(4);
t28 = qJD(1) * t53 + t113;
t24 = qJDD(5) + t28;
t19 = t94 * t24;
t153 = -t104 * t47 + t59 * t19;
t124 = qJD(4) * t67 + t92 * t127 + t95 * t128;
t27 = -qJD(4) * t123 + t124;
t37 = t91 * qJD(4) + t94 * t51;
t9 = t37 * qJD(5) - t94 * qJDD(4) + t91 * t27;
t12 = t95 * t33 - t92 * t34;
t10 = -qJD(4) * pkin(4) - t12;
t110 = t92 * t31 + t95 * t32;
t61 = -t89 * pkin(3) + t70;
t49 = t61 * qJD(1) + qJD(3);
t16 = -t50 * pkin(4) - t51 * pkin(7) + t49;
t120 = qJDD(4) * pkin(7) + t12 * qJD(4) + qJD(5) * t16 + t110;
t148 = pkin(6) + t68;
t54 = t148 * t87;
t55 = t148 * t89;
t107 = -t95 * t54 - t92 * t55;
t14 = -t58 * qJD(3) + t107 * qJD(4);
t18 = t58 * pkin(4) - t59 * pkin(7) + t61;
t21 = -t92 * t54 + t95 * t55;
t152 = -t10 * t52 - (qJD(5) * t18 + t14) * t47 - t120 * t58 + t2 * t59 - t21 * t24;
t151 = g(3) * t76;
t130 = t94 * qJD(4);
t8 = qJD(5) * t130 + t91 * qJDD(4) - t51 * t132 + t94 * t27;
t149 = t8 * t91;
t147 = t37 * t53 + t8 * t58;
t146 = t10 * t59;
t145 = t18 * t24;
t35 = t91 * t51 - t130;
t144 = t35 * t47;
t143 = t37 * t47;
t142 = t37 * t51;
t141 = t51 * t35;
t140 = t77 * t91;
t139 = t77 * t94;
t138 = t79 * t91;
t137 = t79 * t94;
t136 = t91 * t24;
t134 = t87 ^ 2 + t89 ^ 2;
t131 = qJD(5) * t94;
t11 = qJD(4) * pkin(7) + t13;
t46 = t61 * qJDD(1) + qJDD(3);
t6 = t28 * pkin(4) - t27 * pkin(7) + t46;
t119 = qJD(5) * t11 - t6;
t118 = t47 * t94;
t93 = sin(qJ(1));
t96 = cos(qJ(1));
t114 = g(1) * t93 - g(2) * t96;
t112 = -t53 * t35 - t58 * t9;
t38 = -t87 * t56 + t73;
t109 = -t38 * t87 + t39 * t89;
t108 = (-t87 * t62 + t75) * t87 - t45 * t89;
t106 = t19 + (t50 * t91 - t132) * t47;
t105 = -t120 + t151;
t102 = -t129 + t155;
t100 = -pkin(7) * t24 + (t10 + t12) * t47;
t99 = (-t59 * t131 + t91 * t52) * t47 - t59 * t136;
t43 = t78 * t137 + t140;
t42 = -t78 * t138 + t139;
t41 = -t78 * t139 + t138;
t40 = t78 * t140 + t137;
t30 = -t53 * qJD(4) - t58 * qJDD(4);
t29 = -t52 * qJD(4) + t59 * qJDD(4);
t26 = t53 * pkin(4) + t52 * pkin(7);
t15 = t59 * qJD(3) + t21 * qJD(4);
t5 = t94 * t6;
t4 = t94 * t11 + t91 * t16;
t3 = -t91 * t11 + t94 * t16;
t1 = [qJDD(1), t114, g(1) * t96 + g(2) * t93, (t114 + (t88 ^ 2 + t90 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t102 * t89, -t102 * t87, t56 * t134 + t109 - t116, t60 * t70 - g(1) * (-t93 * pkin(1) - t77 * pkin(2) + t79 * qJ(3)) - g(2) * (t96 * pkin(1) + t79 * pkin(2) + t77 * qJ(3)) + t109 * t68 - t108 * qJD(3), t27 * t59 - t51 * t52, -t27 * t58 - t59 * t28 - t52 * t50 - t51 * t53, t29, t30, 0, -t15 * qJD(4) + qJDD(4) * t107 + t115 * t78 + t61 * t28 + t46 * t58 + t49 * t53, -t14 * qJD(4) - t21 * qJDD(4) - t115 * t76 + t61 * t27 + t46 * t59 - t49 * t52, t8 * t94 * t59 - t104 * t37, -(-t35 * t94 - t37 * t91) * t52 + (-t149 - t9 * t94 + (t35 * t91 - t37 * t94) * qJD(5)) * t59, t147 + t153, t112 + t99, t24 * t58 + t47 * t53, -g(1) * t41 - g(2) * t43 + t15 * t35 - t107 * t9 + t3 * t53 + t5 * t58 + (t145 + t26 * t47 + (-t11 * t58 - t21 * t47 + t146) * qJD(5)) * t94 + t152 * t91, -g(1) * t40 - g(2) * t42 + t15 * t37 - t107 * t8 - t4 * t53 + (-(-qJD(5) * t21 + t26) * t47 - t145 + t119 * t58 - qJD(5) * t146) * t91 + t152 * t94; 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, t38 * t89 + t39 * t87 - g(3), 0, 0, 0, 0, 0, t30, -t29, 0, 0, 0, 0, 0, -t112 + t99, t147 - t153; 0, 0, 0, 0, -t127, t128, -t134 * qJD(1) ^ 2, t108 * qJD(1) - t155, 0, 0, 0, 0, 0, 0.2e1 * t51 * qJD(4) + t113, (t50 - t123) * qJD(4) + t124, 0, 0, 0, 0, 0, t106 - t141, -t47 ^ 2 * t94 - t136 - t142; 0, 0, 0, 0, 0, 0, 0, 0, -t51 * t50, -t50 ^ 2 + t51 ^ 2, (-t50 - t123) * qJD(4) + t124, -t113, qJDD(4), -t49 * t51 + t101 - t111, t116 * t78 - t49 * t50 - t110 + t151, t37 * t118 + t149, (t8 - t144) * t94 + (-t9 - t143) * t91, t47 * t118 + t136 - t142, t106 + t141, -t47 * t51, -pkin(4) * t9 + t100 * t91 - t13 * t35 + t156 * t94 - t3 * t51, -pkin(4) * t8 + t100 * t94 - t13 * t37 - t156 * t91 + t4 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37 * t35, -t35 ^ 2 + t37 ^ 2, t8 + t144, t143 - t9, t24, -g(1) * t42 + g(2) * t40 - t10 * t37 + t105 * t91 - t11 * t131 + t4 * t47 + t5, g(1) * t43 - g(2) * t41 + t10 * t35 + t105 * t94 + t119 * t91 + t3 * t47;];
tau_reg = t1;
