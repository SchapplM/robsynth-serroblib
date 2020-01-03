% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RPRP5
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
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPRP5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP5_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:45:04
% EndTime: 2019-12-31 16:45:05
% DurationCPUTime: 0.67s
% Computational Cost: add. (927->170), mult. (2256->204), div. (0->0), fcn. (1604->8), ass. (0->99)
t121 = qJDD(1) * pkin(1);
t76 = sin(qJ(1));
t77 = cos(qJ(1));
t138 = g(1) * t76 - g(2) * t77;
t94 = t138 - qJDD(2) + t121;
t101 = g(1) * t77 + g(2) * t76;
t128 = cos(qJ(3));
t72 = sin(pkin(6));
t73 = cos(pkin(6));
t75 = sin(qJ(3));
t46 = t128 * t72 + t75 * t73;
t137 = t46 * qJD(1);
t132 = t137 ^ 2;
t110 = t128 * t73;
t102 = qJD(1) * t110;
t125 = t75 * t72;
t111 = qJD(1) * t125;
t38 = -t102 + t111;
t35 = t38 ^ 2;
t140 = -t35 - t132;
t139 = -t35 + t132;
t136 = qJ(2) * qJDD(1);
t90 = t110 - t125;
t124 = pkin(5) + qJ(2);
t51 = t124 * t72;
t52 = t124 * t73;
t91 = -t128 * t51 - t75 * t52;
t13 = qJD(2) * t90 + qJD(3) * t91;
t25 = t128 * t52 - t75 * t51;
t71 = pkin(6) + qJ(3);
t66 = sin(t71);
t135 = -t13 * qJD(3) - t25 * qJDD(3) - t138 * t66;
t106 = qJD(3) * t128;
t115 = qJD(1) * qJD(2);
t133 = t124 * qJDD(1) + t115;
t30 = t133 * t72;
t31 = t133 * t73;
t47 = qJD(1) * t51;
t113 = -t47 * t106 + t128 * t31 - t75 * t30;
t67 = cos(t71);
t134 = -g(3) * t66 - t101 * t67 + t113;
t127 = t137 * t38;
t48 = qJD(1) * t52;
t126 = t75 * t48;
t69 = t72 ^ 2;
t70 = t73 ^ 2;
t123 = t69 + t70;
t122 = qJD(3) * t75;
t120 = qJDD(3) * pkin(3);
t23 = t128 * t48 - t75 * t47;
t119 = t23 * qJD(3);
t22 = -t128 * t47 - t126;
t118 = qJD(4) - t22;
t117 = t72 * qJDD(1);
t116 = t73 * qJDD(1);
t114 = qJDD(3) * qJ(4);
t105 = qJDD(1) * t128;
t112 = qJD(3) * t102 + t72 * t105 + t75 * t116;
t64 = t73 * pkin(2) + pkin(1);
t108 = t123 * qJD(1) ^ 2;
t107 = t22 + t126;
t104 = t48 * t106 - t47 * t122 + t128 * t30 + t75 * t31;
t103 = 0.2e1 * t123;
t99 = -t73 * t105 + t117 * t75;
t98 = t67 * pkin(3) + t66 * qJ(4);
t43 = t46 * qJD(3);
t20 = qJD(1) * t43 + t99;
t96 = -t20 * t90 + t38 * t43;
t92 = qJD(3) * t43 - qJDD(3) * t90;
t50 = -qJD(1) * t64 + qJD(2);
t49 = -qJDD(1) * t64 + qJDD(2);
t89 = -g(3) * t67 + t101 * t66 - t104;
t87 = t94 + t121;
t14 = qJD(2) * t46 + qJD(3) * t25;
t86 = -t14 * qJD(3) + t91 * qJDD(3) + t138 * t67;
t19 = qJD(3) * t111 - t112;
t42 = -t106 * t73 + t122 * t72;
t85 = -t137 * t43 - t19 * t90 - t46 * t20 + t42 * t38;
t84 = -t13 * t38 + t137 * t14 + t19 * t91 - t25 * t20 - t101;
t12 = t38 * pkin(3) - qJ(4) * t137 + t50;
t83 = t12 * t137 + qJDD(4) - t89;
t82 = t103 * t115 - t101;
t81 = t20 * pkin(3) + t19 * qJ(4) + t49;
t80 = 0.2e1 * t137 * qJD(3) + t99;
t53 = t77 * t64;
t21 = -t42 * qJD(3) + t46 * qJDD(3);
t18 = -pkin(3) * t90 - t46 * qJ(4) - t64;
t17 = pkin(3) * t137 + t38 * qJ(4);
t16 = qJD(3) * qJ(4) + t23;
t15 = -qJD(3) * pkin(3) + t118;
t11 = t43 * pkin(3) + t42 * qJ(4) - t46 * qJD(4);
t10 = (t38 - t111) * qJD(3) + t112;
t9 = (t38 + t111) * qJD(3) - t112;
t6 = -t137 * t42 - t19 * t46;
t4 = -t122 * t48 + t113;
t3 = qJDD(4) + t104 - t120;
t2 = t114 + (qJD(4) - t126) * qJD(3) + t113;
t1 = -qJD(4) * t137 + t81;
t5 = [0, 0, 0, 0, 0, qJDD(1), t138, t101, 0, 0, t69 * qJDD(1), 0.2e1 * t72 * t116, 0, t70 * qJDD(1), 0, 0, t87 * t73, -t87 * t72, t103 * t136 + t82, t94 * pkin(1) + (t123 * t136 + t82) * qJ(2), t6, t85, t21, t96, -t92, 0, -t64 * t20 + t50 * t43 - t49 * t90 + t86, t64 * t19 - t50 * t42 + t49 * t46 + t135, t104 * t46 + t22 * t42 - t23 * t43 + t4 * t90 + t84, t4 * t25 + t23 * t13 - t104 * t91 - t22 * t14 - t49 * t64 - g(1) * (t124 * t77 - t76 * t64) - g(2) * (t124 * t76 + t53), t6, t21, -t85, 0, t92, t96, -t1 * t90 + t11 * t38 + t12 * t43 + t18 * t20 + t86, -t15 * t42 - t16 * t43 + t2 * t90 + t3 * t46 + t84, -t1 * t46 - t11 * t137 + t12 * t42 + t18 * t19 - t135, -g(2) * t53 + t1 * t18 + t12 * t11 + t16 * t13 + t15 * t14 + t2 * t25 - t3 * t91 + (-g(1) * t124 - g(2) * t98) * t77 + (-g(1) * (-t64 - t98) - g(2) * t124) * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t116, t117, -t108, -qJ(2) * t108 - t94, 0, 0, 0, 0, 0, 0, t80, -t9, t140, t137 * t22 + t23 * t38 - t138 + t49, 0, 0, 0, 0, 0, 0, t80, t140, t9, t16 * t38 + (-qJD(4) - t15) * t137 + t81 - t138; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t127, t139, t10, -t127, -t99, qJDD(3), -t137 * t50 + t119 + t89, qJD(3) * t107 + t50 * t38 - t134, 0, 0, t127, t10, -t139, qJDD(3), t99, -t127, -t17 * t38 + t119 + 0.2e1 * t120 - t83, pkin(3) * t19 - t20 * qJ(4) + (t16 - t23) * t137 + (t15 - t118) * t38, 0.2e1 * t114 - t12 * t38 + t17 * t137 + (0.2e1 * qJD(4) - t107) * qJD(3) + t134, -t3 * pkin(3) - g(3) * t98 + t2 * qJ(4) + t118 * t16 - t12 * t17 - t15 * t23 + t101 * (pkin(3) * t66 - qJ(4) * t67); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(3) + t127, t10, -qJD(3) ^ 2 - t132, -t16 * qJD(3) - t120 + t83;];
tau_reg = t5;
