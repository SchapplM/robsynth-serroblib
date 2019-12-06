% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRPRP1
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
% tau_reg [5x19]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRPRP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:28:42
% EndTime: 2019-12-05 15:28:45
% DurationCPUTime: 0.63s
% Computational Cost: add. (930->168), mult. (1915->198), div. (0->0), fcn. (1389->8), ass. (0->97)
t124 = qJDD(2) * pkin(2);
t79 = pkin(7) + qJ(2);
t73 = sin(t79);
t75 = cos(t79);
t140 = g(1) * t73 - g(2) * t75;
t96 = t140 - qJDD(3) + t124;
t116 = qJD(2) * qJD(3);
t139 = qJ(3) * qJDD(2) + t116;
t104 = g(1) * t75 + g(2) * t73;
t131 = cos(qJ(4));
t80 = sin(pkin(8));
t81 = cos(pkin(8));
t83 = sin(qJ(4));
t45 = t131 * t80 + t83 * t81;
t40 = t45 * qJD(2);
t111 = t131 * t81;
t128 = t83 * t80;
t92 = t111 - t128;
t127 = pkin(6) + qJ(3);
t53 = t127 * t80;
t54 = t127 * t81;
t93 = -t131 * t53 - t83 * t54;
t11 = qJD(3) * t92 + qJD(4) * t93;
t27 = t131 * t54 - t83 * t53;
t78 = pkin(8) + qJ(4);
t72 = sin(t78);
t138 = -t11 * qJD(4) - t27 * qJDD(4) - t140 * t72;
t108 = qJD(4) * t131;
t69 = t81 * qJDD(1);
t28 = t69 + (-t127 * qJDD(2) - t116) * t80;
t118 = t81 * qJDD(2);
t34 = t80 * qJDD(1) + t139 * t81;
t29 = pkin(6) * t118 + t34;
t71 = t81 * qJD(1);
t35 = -qJD(2) * t53 + t71;
t114 = t35 * t108 + t131 * t29 + t83 * t28;
t74 = cos(t78);
t137 = -g(3) * t72 - t104 * t74 + t114;
t136 = t40 ^ 2;
t107 = qJDD(2) * t131;
t119 = t80 * qJDD(2);
t102 = -t81 * t107 + t119 * t83;
t43 = t45 * qJD(4);
t18 = qJD(2) * t43 + t102;
t105 = qJD(2) * t111;
t112 = qJD(2) * t128;
t38 = -t105 + t112;
t125 = qJD(4) * t83;
t42 = -t108 * t81 + t125 * t80;
t132 = -t45 * t18 + t42 * t38;
t130 = t40 * t38;
t120 = qJ(3) * qJD(2);
t47 = t80 * qJD(1) + t81 * t120;
t36 = t81 * qJD(2) * pkin(6) + t47;
t129 = t83 * t36;
t126 = t80 ^ 2 + t81 ^ 2;
t123 = qJDD(4) * pkin(4);
t14 = t131 * t36 + t83 * t35;
t122 = t14 * qJD(4);
t13 = t131 * t35 - t129;
t121 = qJD(5) - t13;
t115 = qJDD(4) * qJ(5);
t113 = qJD(4) * t105 + t80 * t107 + t83 * t118;
t66 = t81 * pkin(3) + pkin(2);
t109 = t13 + t129;
t106 = t36 * t108 + t35 * t125 - t131 * t28 + t83 * t29;
t101 = t74 * pkin(4) + t72 * qJ(5);
t17 = qJD(4) * t112 - t113;
t99 = t17 * t92 + t43 * t40;
t98 = (-t120 * t80 + t71) * t80 - t47 * t81;
t20 = -t43 * qJD(4) + qJDD(4) * t92;
t94 = t101 + t66;
t52 = -qJD(2) * t66 + qJD(3);
t48 = -qJDD(2) * t66 + qJDD(3);
t91 = -g(3) * t74 + t104 * t72 - t106;
t90 = t96 + t124;
t33 = -t139 * t80 + t69;
t89 = -t33 * t80 + t34 * t81 - t104;
t12 = qJD(3) * t45 + qJD(4) * t27;
t88 = -t12 * qJD(4) + t93 * qJDD(4) + t140 * t74;
t9 = t38 * pkin(4) - t40 * qJ(5) + t52;
t87 = t9 * t40 + qJDD(5) - t91;
t86 = t18 * pkin(4) + t17 * qJ(5) + t48;
t85 = 0.2e1 * t40 * qJD(4) + t102;
t37 = t38 ^ 2;
t19 = -t42 * qJD(4) + t45 * qJDD(4);
t16 = -pkin(4) * t92 - t45 * qJ(5) - t66;
t15 = t40 * pkin(4) + t38 * qJ(5);
t10 = qJD(4) * qJ(5) + t14;
t8 = -qJD(4) * pkin(4) + t121;
t6 = t43 * pkin(4) + t42 * qJ(5) - t45 * qJD(5);
t5 = (t38 - t112) * qJD(4) + t113;
t4 = (t38 + t112) * qJD(4) - t113;
t3 = -t40 * qJD(5) + t86;
t2 = qJDD(5) + t106 - t123;
t1 = t115 + (qJD(5) - t129) * qJD(4) + t114;
t7 = [qJDD(1) - g(3), 0, 0, 0, 0, 0, 0, t33 * t81 + t34 * t80 - g(3), 0, 0, 0, 0, 0, t20, -t19, t20, t99 + t132, t19, t1 * t45 - t10 * t42 - t2 * t92 + t8 * t43 - g(3); 0, qJDD(2), t140, t104, t90 * t81, -t90 * t80, t139 * t126 + t89, pkin(2) * t96 + qJ(3) * t89 - qJD(3) * t98, -t17 * t45 - t40 * t42, -t99 + t132, t19, t20, 0, -t66 * t18 + t52 * t43 - t48 * t92 + t88, t66 * t17 - t52 * t42 + t48 * t45 + t138, t16 * t18 - t3 * t92 + t6 * t38 + t9 * t43 + t88, t1 * t92 - t10 * t43 - t11 * t38 + t12 * t40 + t17 * t93 - t27 * t18 + t2 * t45 - t8 * t42 - t104, t16 * t17 - t3 * t45 - t6 * t40 + t9 * t42 - t138, t1 * t27 + t10 * t11 + t8 * t12 + t3 * t16 - t2 * t93 + t9 * t6 + (-g(1) * t127 - g(2) * t94) * t75 + (g(1) * t94 - g(2) * t127) * t73; 0, 0, 0, 0, -t118, t119, -t126 * qJD(2) ^ 2, qJD(2) * t98 - t96, 0, 0, 0, 0, 0, t85, -t4, t85, -t37 - t136, t4, t10 * t38 + (-qJD(5) - t8) * t40 + t86 - t140; 0, 0, 0, 0, 0, 0, 0, 0, t130, -t37 + t136, t5, -t102, qJDD(4), -t52 * t40 + t122 + t91, qJD(4) * t109 + t52 * t38 - t137, -t15 * t38 + t122 + 0.2e1 * t123 - t87, pkin(4) * t17 - t18 * qJ(5) + (t10 - t14) * t40 + (t8 - t121) * t38, 0.2e1 * t115 + t15 * t40 - t9 * t38 + (0.2e1 * qJD(5) - t109) * qJD(4) + t137, -t2 * pkin(4) - g(3) * t101 + t1 * qJ(5) + t10 * t121 - t8 * t14 - t9 * t15 + t104 * (pkin(4) * t72 - qJ(5) * t74); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(4) + t130, t5, -qJD(4) ^ 2 - t136, -t10 * qJD(4) - t123 + t87;];
tau_reg = t7;
