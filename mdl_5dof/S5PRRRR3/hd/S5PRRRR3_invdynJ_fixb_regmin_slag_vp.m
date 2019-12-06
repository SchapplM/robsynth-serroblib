% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRRRR3
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% tau_reg [5x17]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:06
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRRR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:06:21
% EndTime: 2019-12-05 17:06:23
% DurationCPUTime: 0.49s
% Computational Cost: add. (672->121), mult. (962->158), div. (0->0), fcn. (535->12), ass. (0->92)
t50 = pkin(9) + qJ(2);
t48 = qJ(3) + t50;
t38 = qJ(4) + t48;
t33 = sin(t38);
t34 = cos(t38);
t55 = sin(qJ(4));
t56 = sin(qJ(3));
t94 = qJD(4) * t56;
t83 = qJD(2) * t94;
t79 = pkin(2) * t83;
t119 = g(1) * t34 + g(2) * t33 + t55 * t79;
t59 = cos(qJ(3));
t113 = pkin(2) * t59;
t42 = qJDD(2) * t113;
t49 = qJDD(2) + qJDD(3);
t98 = pkin(2) * qJD(2);
t88 = t56 * t98;
t15 = pkin(3) * t49 - qJD(3) * t88 + t42;
t51 = qJD(2) + qJD(3);
t22 = pkin(3) * t51 + t59 * t98;
t58 = cos(qJ(4));
t97 = qJD(3) * t59;
t84 = qJD(2) * t97;
t90 = qJDD(2) * t56;
t66 = (t84 + t90) * pkin(2);
t61 = -(qJD(4) * t22 + t66) * t58 - t55 * t15 + t119;
t93 = qJD(4) * t58;
t87 = t56 * t93;
t122 = ((t55 * t97 + t87) * qJD(2) + t55 * t90) * pkin(2);
t44 = qJDD(4) + t49;
t109 = t44 * pkin(4);
t110 = g(2) * t34;
t13 = t58 * t15;
t95 = qJD(4) * t55;
t76 = t22 * t95 - t13;
t121 = -t109 + t76 + t122 + t110;
t28 = g(1) * t33;
t117 = t28 - t110;
t103 = t56 * t58;
t74 = t55 * t59 + t103;
t18 = t74 * t98;
t39 = pkin(3) * t55 + pkin(8);
t40 = -pkin(3) * t58 - pkin(4);
t47 = qJD(4) + t51;
t60 = qJD(5) ^ 2;
t115 = t39 * t60 + t40 * t44 + (pkin(3) * t95 - t18) * t47;
t112 = pkin(3) * t44;
t111 = pkin(4) * t47;
t41 = pkin(3) + t113;
t108 = (t41 * t95 + (t74 * qJD(3) + t87) * pkin(2)) * t47;
t54 = sin(qJ(5));
t57 = cos(qJ(5));
t11 = t22 * t58 - t55 * t88;
t9 = -t11 - t111;
t107 = t9 * qJD(5) * t54 + t57 * t28;
t106 = (t22 * t55 + t58 * t88) * t47;
t105 = t41 * t58;
t104 = t55 * t56;
t102 = t57 * t44;
t36 = sin(t48);
t37 = cos(t48);
t101 = g(1) * t37 + g(2) * t36;
t100 = pkin(2) * t103 + t55 * t41;
t52 = t54 ^ 2;
t99 = -t57 ^ 2 + t52;
t92 = qJD(5) * t57;
t91 = qJDD(1) - g(3);
t89 = t121 * t54 + t9 * t92;
t81 = qJD(2) * (-qJD(3) + t51);
t80 = qJD(3) * (-qJD(2) - t51);
t78 = g(1) * t36 - g(2) * t37 + t42;
t73 = t58 * t59 - t104;
t71 = pkin(8) * t60 - t106 - t109;
t16 = pkin(2) * t104 - pkin(4) - t105;
t17 = pkin(8) + t100;
t70 = t16 * t44 + t17 * t60 + t108;
t68 = -t44 * pkin(8) - t47 * t9 + t61;
t67 = -t76 + t117;
t65 = -pkin(8) * qJDD(5) + (t11 - t111) * qJD(5);
t4 = t41 * t93 + (t73 * qJD(3) - t55 * t94) * pkin(2);
t64 = -qJDD(5) * t17 + (t16 * t47 - t4) * qJD(5);
t19 = t73 * t98;
t63 = -qJDD(5) * t39 + (-pkin(3) * t93 + t40 * t47 + t19) * qJD(5);
t62 = (-pkin(3) * t47 - t22) * qJD(4) - t66;
t46 = cos(t50);
t45 = sin(t50);
t43 = t47 ^ 2;
t24 = qJDD(5) * t57 - t54 * t60;
t23 = qJDD(5) * t54 + t57 * t60;
t14 = 0.2e1 * t47 * t54 * t92 + t44 * t52;
t8 = -0.2e1 * t99 * t47 * qJD(5) + 0.2e1 * t54 * t102;
t1 = [t91, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t23; 0, qJDD(2), g(1) * t45 - g(2) * t46, g(1) * t46 + g(2) * t45, t49, (t49 * t59 + t56 * t80) * pkin(2) + t78, ((-qJDD(2) - t49) * t56 + t59 * t80) * pkin(2) + t101, t44, t44 * t105 - t108 + (-t58 * t83 + (-t84 + (-qJDD(2) - t44) * t56) * t55) * pkin(2) + t67, -t100 * t44 - t4 * t47 + t61, t14, t8, t23, t24, 0, t64 * t54 + (-t70 - t121) * t57 + t107, t64 * t57 + (t70 - t28) * t54 + t89; 0, 0, 0, 0, t49, t56 * pkin(2) * t81 + t78, (t59 * t81 - t90) * pkin(2) + t101, t44, t18 * t47 + t13 + (-t79 + t112) * t58 + t62 * t55 + t117, t19 * t47 + (-t15 - t112) * t55 + t62 * t58 + t119, t14, t8, t23, t24, 0, t63 * t54 + (-t115 - t121) * t57 + t107, t63 * t57 + (t115 - t28) * t54 + t89; 0, 0, 0, 0, 0, 0, 0, t44, t106 + t67 - t122, t11 * t47 + t61, t14, t8, t23, t24, 0, t65 * t54 + (-t71 - t121) * t57 + t107, t65 * t57 + (t71 - t28) * t54 + t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54 * t43 * t57, t99 * t43, t54 * t44, t102, qJDD(5), t68 * t54 + t91 * t57, -t91 * t54 + t68 * t57;];
tau_reg = t1;
