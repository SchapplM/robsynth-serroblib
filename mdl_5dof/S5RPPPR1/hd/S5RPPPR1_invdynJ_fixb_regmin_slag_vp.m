% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPPPR1
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% tau_reg [5x19]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPPR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:20:36
% EndTime: 2020-01-03 11:20:41
% DurationCPUTime: 0.91s
% Computational Cost: add. (865->180), mult. (1783->270), div. (0->0), fcn. (1301->14), ass. (0->122)
t81 = qJ(1) + pkin(7);
t71 = sin(t81);
t73 = cos(t81);
t154 = g(2) * t73 + g(3) * t71;
t87 = cos(pkin(7));
t69 = -t87 * pkin(1) - pkin(2);
t129 = qJDD(1) * t69;
t56 = qJDD(3) + t129;
t155 = t154 + t56;
t82 = sin(pkin(9));
t85 = cos(pkin(9));
t88 = sin(qJ(5));
t90 = cos(qJ(5));
t102 = t90 * t82 + t88 * t85;
t83 = sin(pkin(8));
t135 = qJD(1) * t83;
t121 = t82 * t135;
t116 = t88 * t121;
t140 = t90 * t85;
t120 = qJD(5) * t140;
t11 = -qJD(5) * t116 + (qJD(1) * t120 + t102 * qJDD(1)) * t83;
t86 = cos(pkin(8));
t62 = t86 * qJD(1) - qJD(5);
t153 = qJD(5) + t62;
t84 = sin(pkin(7));
t64 = t84 * pkin(1) + qJ(3);
t52 = qJD(1) * qJD(3) + qJDD(1) * t64;
t151 = pkin(6) * t83;
t150 = g(1) * t83;
t123 = t83 * t140;
t128 = qJDD(1) * t82;
t96 = t102 * qJD(5);
t10 = qJDD(1) * t123 + (-qJD(1) * t96 - t88 * t128) * t83;
t149 = t10 * t86;
t35 = t102 * t135;
t148 = t35 * t62;
t37 = qJD(1) * t123 - t116;
t147 = t37 * t62;
t146 = t64 * t86;
t145 = t71 * t86;
t144 = t73 * t86;
t143 = t86 * t11;
t92 = qJD(1) ^ 2;
t142 = t86 * t92;
t141 = t88 * t82;
t39 = (-qJD(5) * t141 + t120) * t83;
t43 = t102 * t83;
t126 = t86 * qJDD(1);
t61 = -qJDD(5) + t126;
t139 = t39 * t62 + t43 * t61;
t132 = qJD(4) * t83;
t51 = -t86 * pkin(3) - t83 * qJ(4) + t69;
t19 = -qJD(1) * t132 + t51 * qJDD(1) + qJDD(3);
t32 = t83 * qJDD(2) + t86 * t52;
t6 = t82 * t19 + t85 * t32;
t33 = t51 * qJD(1) + qJD(3);
t60 = t64 * qJD(1);
t42 = t83 * qJD(2) + t86 * t60;
t9 = t82 * t33 + t85 * t42;
t18 = t85 * t146 + t82 * t51;
t137 = -t82 ^ 2 - t85 ^ 2;
t77 = t83 ^ 2;
t136 = t86 ^ 2 + t77;
t134 = qJD(3) * t83;
t133 = qJD(3) * t86;
t131 = qJDD(2) - g(1);
t127 = t83 * qJDD(1);
t124 = t85 * t151;
t91 = cos(qJ(1));
t122 = t91 * pkin(1) + t73 * pkin(2) + t71 * qJ(3);
t100 = -t86 * pkin(4) - t124;
t5 = t85 * t19 - t82 * t32;
t2 = t100 * qJDD(1) + t5;
t118 = t82 * t127;
t3 = -pkin(6) * t118 + t6;
t119 = t90 * t2 - t88 * t3;
t117 = t136 * t92;
t8 = t85 * t33 - t82 * t42;
t41 = t86 * qJD(2) - t83 * t60;
t48 = t83 * t52;
t31 = t86 * qJDD(2) - t48;
t89 = sin(qJ(1));
t114 = t89 * pkin(1) + t71 * pkin(2) - t73 * qJ(3);
t112 = -g(2) * t71 + g(3) * t73;
t111 = -g(2) * t91 - g(3) * t89;
t110 = t88 * t2 + t90 * t3;
t4 = t100 * qJD(1) + t8;
t7 = -pkin(6) * t121 + t9;
t109 = t90 * t4 - t88 * t7;
t108 = -t88 * t4 - t90 * t7;
t46 = t85 * t51;
t12 = -t124 + t46 + (-t64 * t82 - pkin(4)) * t86;
t13 = -t82 * t151 + t18;
t107 = t90 * t12 - t88 * t13;
t106 = t88 * t12 + t90 * t13;
t105 = -t31 * t83 + t32 * t86;
t38 = t83 * t96;
t101 = -t140 + t141;
t44 = t101 * t83;
t104 = t38 * t62 + t44 * t61;
t103 = t41 * t83 - t42 * t86;
t40 = qJD(4) - t41;
t26 = qJDD(4) - t31;
t99 = t101 * t62;
t17 = -t82 * t146 + t46;
t49 = -t85 * t132 - t82 * t133;
t98 = -t49 * qJD(1) - t17 * qJDD(1) - t5;
t50 = -t82 * t132 + t85 * t133;
t97 = t50 * qJD(1) + t18 * qJDD(1) + t6;
t95 = t129 + t155;
t94 = t26 * t83 + t52 * t77 + t112;
t80 = pkin(9) + qJ(5);
t72 = cos(t80);
t70 = sin(t80);
t47 = (pkin(4) * t82 + t64) * t83;
t30 = t72 * t144 + t71 * t70;
t29 = t70 * t144 - t71 * t72;
t28 = t72 * t145 - t73 * t70;
t27 = -t70 * t145 - t73 * t72;
t21 = pkin(4) * t121 + t40;
t16 = pkin(4) * t118 + t26;
t1 = [qJDD(1), t111, g(2) * t89 - g(3) * t91, (t111 + (t84 ^ 2 + t87 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), -t95 * t86, t95 * t83, t52 * t136 + t105 + t112, -g(2) * t122 - g(3) * t114 - t103 * qJD(3) + t105 * t64 + t56 * t69, (-t154 * t85 + t98) * t86 + t94 * t82, (t154 * t82 + t97) * t86 + t94 * t85, (-t97 * t82 + t98 * t85 - t154) * t83, t6 * t18 + t9 * t50 + t5 * t17 + t8 * t49 - g(2) * (pkin(3) * t144 + t122) - g(3) * (pkin(3) * t145 + t114) + (-qJ(4) * t154 + t40 * qJD(3) + t26 * t64) * t83, -t10 * t44 - t37 * t38, -t10 * t43 + t44 * t11 + t38 * t35 - t37 * t39, t104 - t149, t139 + t143, t61 * t86, -(t90 * t49 - t88 * t50) * t62 - t107 * t61 - t119 * t86 + t35 * t134 + t47 * t11 + t16 * t43 + t21 * t39 - g(2) * t30 - g(3) * t28 + (t106 * t62 - t108 * t86) * qJD(5), (t88 * t49 + t90 * t50) * t62 + t106 * t61 + t110 * t86 + t37 * t134 + t47 * t10 - t16 * t44 - t21 * t38 + g(2) * t29 - g(3) * t27 + (t107 * t62 + t109 * t86) * qJD(5); 0, 0, 0, t131, 0, 0, 0, t31 * t86 + t32 * t83 - g(1), 0, 0, 0, -t26 * t86 - g(1) + (-t5 * t82 + t6 * t85) * t83, 0, 0, 0, 0, 0, t139 - t143, -t104 - t149; 0, 0, 0, 0, -t126, t127, -t117, t103 * qJD(1) + t155, -t82 * t117 - t85 * t126, -t85 * t117 + t82 * t126, t137 * t127, t5 * t85 + t6 * t82 + (-t40 * t83 + (t8 * t82 - t85 * t9) * t86) * qJD(1) + t154, 0, 0, 0, 0, 0, t101 * t61 + t62 * t96 + (-t102 * t62 * t86 - t83 * t35) * qJD(1), t102 * t61 - qJD(5) * t99 + (-t83 * t37 + t86 * t99) * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, (-t85 * t142 + t128) * t83, (qJDD(1) * t85 + t82 * t142) * t83, t137 * t92 * t77, qJDD(4) + t48 - t131 * t86 + ((t8 * t85 + t82 * t9) * qJD(1) + t112) * t83, 0, 0, 0, 0, 0, t11 - t147, t10 + t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37 * t35, -t35 ^ 2 + t37 ^ 2, t10 - t148, -t11 - t147, -t61, -g(2) * t27 - g(3) * t29 + t153 * t108 + t70 * t150 - t21 * t37 + t119, g(2) * t28 - g(3) * t30 - t153 * t109 + t72 * t150 + t21 * t35 - t110;];
tau_reg = t1;
