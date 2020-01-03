% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPPR3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% tau_reg [5x17]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPPR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:26:37
% EndTime: 2019-12-31 19:26:38
% DurationCPUTime: 0.43s
% Computational Cost: add. (619->135), mult. (948->161), div. (0->0), fcn. (531->12), ass. (0->96)
t67 = cos(qJ(2));
t100 = pkin(1) * qJD(2);
t88 = qJD(1) * t100;
t64 = sin(qJ(2));
t97 = qJDD(1) * t64;
t118 = pkin(1) * t97 + t67 * t88;
t56 = qJD(1) + qJD(2);
t101 = pkin(1) * qJD(1);
t90 = t67 * t101;
t25 = pkin(2) * t56 + t90;
t61 = sin(pkin(8));
t62 = cos(pkin(8));
t91 = t64 * t101;
t14 = t25 * t61 + t62 * t91;
t11 = qJ(4) * t56 + t14;
t60 = qJ(1) + qJ(2);
t50 = pkin(8) + t60;
t41 = sin(t50);
t42 = cos(t50);
t112 = pkin(1) * t67;
t49 = qJDD(1) * t112;
t55 = qJDD(1) + qJDD(2);
t18 = pkin(2) * t55 - t64 * t88 + t49;
t6 = -t118 * t61 + t62 * t18;
t76 = qJDD(4) - t6;
t72 = -g(1) * t41 + g(2) * t42 + t76;
t117 = -t11 * t56 + t72;
t51 = sin(t60);
t52 = cos(t60);
t116 = g(1) * t51 - g(2) * t52;
t54 = t56 ^ 2;
t7 = t118 * t62 + t61 * t18;
t4 = t55 * qJ(4) + t56 * qJD(4) + t7;
t63 = sin(qJ(5));
t66 = cos(qJ(5));
t99 = qJD(5) * t66;
t115 = t11 * t99 + t4 * t63;
t114 = pkin(1) * t64;
t65 = sin(qJ(1));
t113 = pkin(1) * t65;
t111 = pkin(2) * t51;
t110 = pkin(3) * t55;
t106 = t62 * t64;
t105 = t66 * t55;
t104 = g(1) * t52 + g(2) * t51;
t69 = qJD(5) ^ 2;
t103 = -t54 - t69;
t58 = t66 ^ 2;
t102 = t63 ^ 2 - t58;
t32 = t61 * t91;
t24 = t62 * t90 - t32;
t98 = qJD(4) - t24;
t59 = qJDD(3) - g(3);
t38 = t61 * t114;
t48 = pkin(2) + t112;
t85 = t48 * t62 - t38;
t21 = -pkin(3) - t85;
t17 = -pkin(7) + t21;
t96 = qJDD(5) * t17;
t43 = -pkin(2) * t62 - pkin(3);
t39 = -pkin(7) + t43;
t95 = qJDD(5) * t39;
t94 = qJDD(5) * t63;
t93 = qJDD(5) * t66;
t47 = pkin(2) * t52;
t92 = t42 * pkin(3) + t41 * qJ(4) + t47;
t77 = pkin(1) * t106 + t48 * t61;
t20 = qJ(4) + t77;
t75 = pkin(1) * (t61 * t67 + t106);
t23 = qJD(2) * t75;
t87 = t20 * t56 + t23;
t22 = qJD(1) * t75;
t40 = pkin(2) * t61 + qJ(4);
t86 = t40 * t56 - t22;
t13 = t62 * t25 - t32;
t84 = qJD(1) * (-qJD(2) + t56);
t83 = qJD(2) * (-qJD(1) - t56);
t81 = t49 + t116;
t80 = -g(1) * t42 - g(2) * t41;
t79 = t62 * t67 * t100 - qJD(2) * t38;
t78 = -pkin(3) * t41 + t42 * qJ(4) - t111;
t74 = t80 + t7;
t73 = -(-pkin(3) - pkin(7)) * t55 - t117;
t16 = qJD(4) + t79;
t71 = t16 * t56 - t17 * t69 + t20 * t55 + t80;
t70 = -t39 * t69 + t40 * t55 + t98 * t56 + t80;
t68 = cos(qJ(1));
t53 = t68 * pkin(1);
t27 = -t63 * t69 + t93;
t26 = -t66 * t69 - t94;
t19 = -0.2e1 * t56 * t63 * t99 + t55 * t58;
t12 = 0.2e1 * t102 * t56 * qJD(5) - 0.2e1 * t63 * t105;
t10 = -pkin(3) * t56 + qJD(4) - t13;
t5 = t76 - t110;
t2 = t4 * t66;
t1 = [qJDD(1), g(1) * t65 - g(2) * t68, g(1) * t68 + g(2) * t65, t55, (t55 * t67 + t64 * t83) * pkin(1) + t81, ((-qJDD(1) - t55) * t64 + t67 * t83) * pkin(1) + t104, t7 * t77 + t14 * t79 + t6 * t85 - t13 * t23 - g(1) * (-t111 - t113) - g(2) * (t47 + t53), t23 * t56 + (-pkin(3) + t21) * t55 + t72, (qJD(4) + t16) * t56 + (qJ(4) + t20) * t55 + t74, t4 * t20 + t11 * t16 + t5 * t21 + t10 * t23 - g(1) * (t78 - t113) - g(2) * (t53 + t92), t19, t12, t27, t26, 0, (t87 * qJD(5) + t96) * t66 + t71 * t63 + t115, t2 + (-t96 + (-t11 - t87) * qJD(5)) * t63 + t71 * t66; 0, 0, 0, t55, t84 * t114 + t81, (t67 * t84 - t97) * pkin(1) + t104, t13 * t22 - t14 * t24 + (t6 * t62 + t61 * t7 + t116) * pkin(2), -t22 * t56 + (-pkin(3) + t43) * t55 + t72, (0.2e1 * qJD(4) - t24) * t56 + (qJ(4) + t40) * t55 + t74, -g(1) * t78 - g(2) * t92 - t10 * t22 + t98 * t11 + t4 * t40 + t5 * t43, t19, t12, t27, t26, 0, (t86 * qJD(5) + t95) * t66 + t70 * t63 + t115, t2 + (-t95 + (-t11 - t86) * qJD(5)) * t63 + t70 * t66; 0, 0, 0, 0, 0, 0, t59, 0, 0, t59, 0, 0, 0, 0, 0, t26, -t27; 0, 0, 0, 0, 0, 0, 0, t55, -t54, -t110 + t117, 0, 0, 0, 0, 0, t103 * t63 + t93, t103 * t66 - t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66 * t54 * t63, -t102 * t54, t105, -t63 * t55, qJDD(5), -t59 * t63 - t73 * t66, -t59 * t66 + t73 * t63;];
tau_reg = t1;
