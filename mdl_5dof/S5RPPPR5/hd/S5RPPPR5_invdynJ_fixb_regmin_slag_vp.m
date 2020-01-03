% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPPPR5
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
% 
% Output:
% tau_reg [5x20]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPPR5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR5_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:46:32
% EndTime: 2019-12-31 17:46:34
% DurationCPUTime: 0.58s
% Computational Cost: add. (622->155), mult. (1106->196), div. (0->0), fcn. (766->10), ass. (0->94)
t75 = sin(pkin(8));
t77 = cos(pkin(8));
t80 = sin(qJ(5));
t81 = cos(qJ(5));
t32 = t81 * t75 + t80 * t77;
t130 = t32 * qJD(1);
t131 = t130 * qJD(5);
t119 = t75 ^ 2 + t77 ^ 2;
t123 = t80 * t75;
t30 = -t81 * t77 + t123;
t28 = t30 * qJD(5);
t129 = qJDD(1) * pkin(3) + qJDD(4);
t128 = qJD(5) ^ 2;
t82 = -pkin(1) - pkin(2);
t76 = sin(pkin(7));
t78 = cos(pkin(7));
t38 = t78 * qJ(2) + t76 * t82;
t34 = -qJ(4) + t38;
t127 = pkin(6) - t34;
t110 = qJ(2) * qJDD(1);
t41 = t82 * qJDD(1) + qJDD(2);
t122 = t78 * t110 + t76 * t41;
t109 = qJD(1) * qJD(2);
t49 = t78 * t109;
t17 = t49 + t122;
t11 = -qJDD(1) * qJ(4) - qJD(1) * qJD(4) + t17;
t4 = t75 * qJDD(3) + t77 * t11;
t126 = cos(qJ(1));
t125 = sin(qJ(1));
t83 = qJD(1) ^ 2;
t124 = t78 * t83;
t115 = qJ(2) * qJD(1);
t42 = t82 * qJD(1) + qJD(2);
t24 = t78 * t115 + t76 * t42;
t121 = t126 * pkin(1) + t125 * qJ(2);
t120 = g(1) * t125 - g(2) * t126;
t118 = pkin(1) * qJDD(1);
t117 = qJD(1) * t77;
t116 = qJD(2) * t76;
t114 = t30 * qJDD(5);
t113 = t32 * qJDD(5);
t112 = t75 * qJDD(1);
t111 = t77 * qJDD(1);
t108 = t126 * pkin(2) + t121;
t107 = qJD(1) * t123;
t106 = t81 * t117;
t105 = 0.2e1 * t109;
t47 = t76 * t109;
t104 = -t76 * t110 + t78 * t41;
t23 = -t76 * t115 + t78 * t42;
t37 = -t76 * qJ(2) + t78 * t82;
t103 = qJDD(1) * t119;
t102 = qJDD(2) - t118;
t35 = pkin(3) - t37;
t16 = t104 - t47;
t101 = -t125 * pkin(1) + t126 * qJ(2);
t31 = -t125 * t76 - t126 * t78;
t33 = -t125 * t78 + t126 * t76;
t100 = g(1) * t33 - g(2) * t31;
t99 = -g(1) * t31 - g(2) * t33;
t58 = t77 * qJDD(3);
t3 = -t75 * t11 + t58;
t98 = -t3 * t75 + t4 * t77;
t97 = -t81 * t111 + t112 * t80;
t96 = t119 * (-qJD(1) * qJ(4) + t24);
t21 = t127 * t75;
t22 = t127 * t77;
t95 = t81 * t21 + t80 * t22;
t94 = t80 * t21 - t81 * t22;
t93 = t23 * t76 - t24 * t78;
t92 = t78 * qJDD(1) + t76 * t83;
t19 = qJD(1) * pkin(3) + qJD(4) - t23;
t91 = t28 * qJD(5) - t113;
t29 = t32 * qJD(5);
t90 = t29 * qJD(5) + t114;
t14 = -t16 + t129;
t88 = g(1) * t126 + g(2) * t125;
t87 = -t100 - t104;
t86 = -t125 * pkin(2) + t101;
t85 = qJD(5) * t107 - qJDD(1) * t32;
t84 = qJDD(1) * t35 - t100 + t14 + t47;
t73 = pkin(8) + qJ(5);
t62 = cos(t73);
t61 = sin(t73);
t45 = t78 * qJD(2) - qJD(4);
t26 = -t106 + t107;
t25 = t77 * pkin(4) + t35;
t15 = pkin(4) * t117 + t19;
t10 = qJD(1) * t29 + t97;
t9 = -qJD(5) * t106 + t85;
t7 = pkin(4) * t111 + t14;
t2 = -pkin(6) * t111 + t4;
t1 = t58 + (pkin(6) * qJDD(1) - t11) * t75;
t5 = [qJDD(1), t120, t88, -qJDD(2) + 0.2e1 * t118 + t120, t105 - t88 + 0.2e1 * t110, -t102 * pkin(1) - g(1) * t101 - g(2) * t121 + (t105 + t110) * qJ(2), -t37 * qJDD(1) + 0.2e1 * t47 + t87, t38 * qJDD(1) + t122 + 0.2e1 * t49 - t99, -g(1) * t86 - g(2) * t108 - qJD(2) * t93 + t16 * t37 + t17 * t38, t84 * t77, -t84 * t75, -qJD(1) * t119 * t45 - t103 * t34 - t98 + t99, t14 * t35 + t19 * t116 - g(1) * (t33 * pkin(3) + t31 * qJ(4) + t86) - g(2) * (-t31 * pkin(3) + t33 * qJ(4) + t108) + t96 * t45 + t98 * t34, -t130 * t28 - t9 * t32, -t32 * t10 - t130 * t29 + t28 * t26 + t9 * t30, t91, t90, 0, -t26 * t116 - t25 * t10 - t7 * t30 - t15 * t29 + t95 * qJDD(5) - t100 * t62 + (-qJD(5) * t94 - t32 * t45) * qJD(5), -t130 * t116 + t25 * t9 - t7 * t32 + t15 * t28 - t94 * qJDD(5) + t100 * t61 + (-qJD(5) * t95 + t30 * t45) * qJD(5); 0, 0, 0, -qJDD(1), -t83, -t83 * qJ(2) + t102 - t120, -t92, t76 * qJDD(1) - t124, t93 * qJD(1) + t16 * t78 + t17 * t76 - t120, -t92 * t77, t92 * t75, -t103 * t76 + t119 * t124, -t14 * t78 + t98 * t76 + (-t19 * t76 - t78 * t96) * qJD(1) - t120, 0, 0, 0, 0, 0, (t10 + t131) * t78 + (qJD(1) * t26 + t128 * t30 - t113) * t76, (-qJD(1) * t28 - t9) * t78 + (qJD(1) * t130 + t128 * t32 + t114) * t76; 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3) + g(3), 0, 0, 0, t3 * t77 + t4 * t75 + g(3), 0, 0, 0, 0, 0, -t90, t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, t111, -t112, -t119 * t83, qJD(1) * t96 + t129 + t47 + t87, 0, 0, 0, 0, 0, -t97 - 0.2e1 * t131, (t26 - t106) * qJD(5) + t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130 * t26, t130 ^ 2 - t26 ^ 2, (-t26 - t106) * qJD(5) + t85, t97, qJDD(5), g(3) * t62 + t81 * t1 + t130 * t15 - t80 * t2 + t61 * t99, -g(3) * t61 - t80 * t1 - t15 * t26 - t81 * t2 + t62 * t99;];
tau_reg = t5;
