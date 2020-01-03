% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPPRP4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% 
% Output:
% tau_reg [5x18]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRP4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP4_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:52:22
% EndTime: 2019-12-31 17:52:24
% DurationCPUTime: 0.65s
% Computational Cost: add. (666->159), mult. (1132->186), div. (0->0), fcn. (632->6), ass. (0->95)
t63 = -qJ(5) - pkin(6);
t66 = -pkin(1) - pkin(2);
t31 = t66 * qJD(1) + qJD(2);
t60 = sin(pkin(7));
t61 = cos(pkin(7));
t99 = qJ(2) * qJD(1);
t19 = t60 * t31 + t61 * t99;
t121 = t63 * qJD(1) + t19;
t68 = qJD(1) ^ 2;
t120 = t60 * qJDD(1) - t61 * t68;
t65 = cos(qJ(4));
t93 = t65 * qJDD(1);
t119 = pkin(4) * t93 + qJDD(5);
t106 = qJD(4) * pkin(4);
t64 = sin(qJ(4));
t7 = t65 * qJD(3) - t121 * t64;
t4 = t7 + t106;
t118 = t4 - t7;
t117 = g(3) * t65;
t8 = t64 * qJD(3) + t121 * t65;
t116 = t65 * t8;
t115 = cos(qJ(1));
t114 = sin(qJ(1));
t30 = t66 * qJDD(1) + qJDD(2);
t92 = qJ(2) * qJDD(1);
t112 = t60 * t30 + t61 * t92;
t26 = t61 * qJ(2) + t60 * t66;
t111 = t115 * pkin(1) + t114 * qJ(2);
t110 = g(1) * t114 - g(2) * t115;
t57 = t64 ^ 2;
t58 = t65 ^ 2;
t109 = t57 - t58;
t108 = t57 + t58;
t67 = qJD(4) ^ 2;
t107 = t67 + t68;
t23 = -pkin(6) + t26;
t105 = qJ(5) - t23;
t104 = pkin(1) * qJDD(1);
t102 = qJD(1) * t65;
t18 = t61 * t31 - t60 * t99;
t14 = qJD(1) * pkin(3) - t18;
t11 = pkin(4) * t102 + qJD(5) + t14;
t103 = qJD(1) * t11;
t101 = qJD(2) * t61;
t100 = qJDD(3) + g(3);
t97 = qJDD(4) * t64;
t96 = qJDD(4) * t65;
t94 = t64 * qJDD(1);
t91 = qJD(1) * qJD(2);
t90 = qJD(1) * qJD(4);
t35 = t61 * t91;
t13 = t35 + t112;
t89 = t115 * pkin(2) + t111;
t88 = 0.2e1 * t91;
t87 = t64 * t90;
t33 = t60 * t91;
t86 = t61 * t30 - t60 * t92;
t25 = -t60 * qJ(2) + t61 * t66;
t85 = qJD(4) * t105;
t84 = 0.2e1 * t65 * t90;
t83 = -qJD(5) + t101;
t82 = qJDD(2) - t104;
t22 = pkin(3) - t25;
t12 = -t33 + t86;
t81 = -t114 * pkin(1) + t115 * qJ(2);
t20 = -t114 * t60 - t115 * t61;
t21 = -t114 * t61 + t115 * t60;
t80 = -g(1) * t21 + g(2) * t20;
t79 = -g(1) * t20 - g(2) * t21;
t78 = t4 * t64 - t116;
t77 = t18 * t60 - t19 * t61;
t59 = qJDD(1) * pkin(3);
t9 = -t12 + t59;
t76 = t121 * qJD(4);
t75 = g(1) * t115 + g(2) * t114;
t74 = t80 - t86;
t73 = -t114 * pkin(2) + t81;
t10 = -qJDD(1) * pkin(6) + t13;
t72 = t14 * qJD(1) - t10 + t79;
t71 = qJ(5) * qJDD(1) + qJD(1) * qJD(5) - qJD(4) * qJD(3) - t10;
t70 = -qJDD(4) * t23 + (-qJD(1) * t22 - t101 - t14) * qJD(4);
t69 = qJDD(1) * t22 - t23 * t67 + t33 + t80 + t9;
t52 = t65 * pkin(4);
t48 = t65 * qJDD(3);
t44 = t52 + pkin(3);
t29 = -t67 * t64 + t96;
t28 = -t67 * t65 - t97;
t17 = t105 * t65;
t16 = t105 * t64;
t6 = -t83 * t64 + t65 * t85;
t5 = t64 * t85 + t83 * t65;
t3 = -pkin(4) * t87 + t119 + t9;
t2 = (qJDD(3) - t76) * t64 - t71 * t65;
t1 = qJDD(4) * pkin(4) + t71 * t64 - t65 * t76 + t48;
t15 = [qJDD(1), t110, t75, -qJDD(2) + 0.2e1 * t104 + t110, -t75 + t88 + 0.2e1 * t92, -t82 * pkin(1) - g(1) * t81 - g(2) * t111 + (t88 + t92) * qJ(2), -t25 * qJDD(1) + 0.2e1 * t33 + t74, t26 * qJDD(1) + t112 + 0.2e1 * t35 - t79, -g(1) * t73 - g(2) * t89 - t77 * qJD(2) + t12 * t25 + t13 * t26, t57 * qJDD(1) + t64 * t84, -0.2e1 * t109 * t90 + 0.2e1 * t64 * t93, t28, -t29, 0, t70 * t64 + t69 * t65, -t69 * t64 + t70 * t65, (qJD(4) * t4 + qJDD(1) * t17 - t2 + (qJD(4) * t16 - t5) * qJD(1)) * t65 + (qJD(4) * t8 + qJDD(1) * t16 + t1 + (-qJD(4) * t17 + t6) * qJD(1)) * t64 + t79, -t2 * t17 + t8 * t5 + t1 * t16 + t4 * t6 + t3 * (t22 + t52) + t11 * (t60 * qJD(2) - t64 * t106) - g(1) * (-t20 * t63 + t21 * t44 + t73) - g(2) * (-t20 * t44 - t21 * t63 + t89); 0, 0, 0, -qJDD(1), -t68, -t68 * qJ(2) - t110 + t82, -t61 * qJDD(1) - t60 * t68, t120, t77 * qJD(1) + t12 * t61 + t13 * t60 - t110, 0, 0, 0, 0, 0, (0.2e1 * t87 - t93) * t61 + (-t107 * t65 - t97) * t60, (t84 + t94) * t61 + (t107 * t64 - t96) * t60, -t120 * t108, (t78 * qJD(1) - t3) * t61 + (-t103 - t1 * t64 + t2 * t65 + (-t4 * t65 - t64 * t8) * qJD(4)) * t60 - t110; 0, 0, 0, 0, 0, 0, 0, 0, t100, 0, 0, 0, 0, 0, t29, t28, 0, -t78 * qJD(4) + t1 * t65 + t2 * t64 + g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64 * t68 * t65, t109 * t68, -t94, -t93, qJDD(4), t72 * t64 + t117 + t48, -t100 * t64 + t72 * t65, pkin(4) * t94 + (t106 - t118) * t102, t118 * t8 + (t117 + t1 + (t79 + t103) * t64) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108 * t68, t33 + t59 + (t116 + (-t4 - t106) * t64) * qJD(1) + t74 + t119;];
tau_reg = t15;
