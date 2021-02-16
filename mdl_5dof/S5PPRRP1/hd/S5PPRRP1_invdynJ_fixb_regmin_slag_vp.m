% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PPRRP1
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
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% tau_reg [5x16]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:48
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PPRRP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:48:17
% EndTime: 2021-01-15 14:48:21
% DurationCPUTime: 0.78s
% Computational Cost: add. (670->161), mult. (1467->202), div. (0->0), fcn. (1080->10), ass. (0->99)
t54 = sin(pkin(8));
t56 = cos(pkin(8));
t60 = sin(qJ(3));
t62 = cos(qJ(3));
t131 = -t60 * t54 + t62 * t56;
t24 = t131 * qJD(1);
t51 = pkin(8) + qJ(3);
t46 = cos(t51);
t122 = g(3) * t46;
t45 = sin(t51);
t55 = sin(pkin(7));
t57 = cos(pkin(7));
t82 = g(1) * t57 + g(2) * t55;
t76 = t82 * t45;
t132 = t76 - t122;
t29 = t62 * t54 + t60 * t56;
t27 = t29 * qJD(3);
t61 = cos(qJ(4));
t25 = t29 * qJD(1);
t17 = qJD(3) * pkin(6) + t25;
t96 = qJ(5) * qJD(3);
t85 = t17 + t96;
t79 = t85 * t61;
t103 = qJD(3) * t131;
t130 = qJD(3) * t103 + t29 * qJDD(3);
t90 = -g(1) * t55 + g(2) * t57;
t129 = t29 * qJDD(1);
t123 = g(3) * t45;
t70 = t82 * t46 + t123;
t104 = qJD(4) * pkin(4);
t49 = t61 * qJD(2);
t59 = sin(qJ(4));
t9 = -t85 * t59 + t49;
t8 = t9 + t104;
t128 = -t9 + t8;
t52 = t59 ^ 2;
t127 = pkin(4) * t52;
t121 = g(3) * t59;
t120 = t61 * pkin(4);
t119 = t55 * t59;
t118 = t55 * t61;
t117 = t57 * t59;
t116 = t57 * t61;
t113 = t61 * t24;
t64 = qJD(3) ^ 2;
t112 = t61 * t64;
t109 = qJ(5) + pkin(6);
t108 = qJD(4) * t113 + t46 * t121;
t53 = t61 ^ 2;
t107 = t52 - t53;
t106 = t52 + t53;
t105 = qJD(3) * pkin(3);
t101 = qJDD(4) * pkin(4);
t16 = -t24 - t105;
t100 = t16 * qJD(3);
t99 = t25 * qJD(3);
t97 = t59 * qJD(4);
t44 = pkin(3) + t120;
t95 = qJDD(3) * t44;
t47 = t59 * qJDD(3);
t93 = t61 * qJDD(3);
t91 = qJD(3) * qJD(4);
t89 = t59 * t91;
t88 = qJD(4) * t109;
t87 = t24 * t97 + t61 * t99 + (g(1) * t116 + g(2) * t118) * t45;
t74 = -qJD(1) * t27 + t131 * qJDD(1);
t5 = pkin(4) * t89 + qJDD(5) - t74 - t95;
t86 = t5 - t95;
t6 = qJDD(3) * pkin(6) + qJD(1) * t103 + t129;
t84 = -qJD(4) * qJD(2) - t6;
t83 = 0.2e1 * t61 * t91;
t10 = t59 * qJD(2) + t79;
t81 = t10 * t61 - t8 * t59;
t78 = -t27 * qJD(3) + qJDD(3) * t131;
t63 = qJD(4) ^ 2;
t77 = -0.2e1 * qJDD(3) * pkin(3) + pkin(6) * t63 - t74;
t48 = t61 * qJDD(2);
t75 = -g(1) * (-t46 * t117 + t118) - g(2) * (-t46 * t119 - t116) + t45 * t121 + t48;
t73 = -qJ(5) * qJDD(3) + t84;
t72 = t29 * t63 - t78;
t71 = -pkin(6) * qJDD(4) + (t16 - t105) * qJD(4);
t69 = -0.2e1 * t103 * qJD(4) - qJDD(4) * t29;
t12 = t17 * t97;
t68 = -g(1) * (-t46 * t116 - t119) - g(2) * (-t46 * t118 + t117) - t59 * qJDD(2) + t12 + t61 * t123;
t67 = -t76 - t99;
t66 = qJD(3) * qJD(5) - t73;
t11 = -t44 * qJD(3) + qJD(5) - t24;
t65 = (-qJD(5) - t11) * qJD(3) + t73;
t35 = t109 * t61;
t34 = t109 * t59;
t33 = qJDD(4) * t61 - t63 * t59;
t32 = qJDD(4) * t59 + t63 * t61;
t23 = -t59 * qJD(5) - t61 * t88;
t22 = t61 * qJD(5) - t59 * t88;
t4 = -t12 + (-qJ(5) * t91 + qJDD(2)) * t59 + t66 * t61;
t3 = -qJD(4) * t79 - t66 * t59 + t101 + t48;
t2 = t69 * t59 - t72 * t61;
t1 = t72 * t59 + t69 * t61;
t7 = [qJDD(1) - g(3), -g(3) + (t54 ^ 2 + t56 ^ 2) * qJDD(1), 0, t78, -t130, 0, 0, 0, 0, 0, t2, t1, t2, t1, t130 * t106, t11 * t27 - t5 * t131 - g(3) + t81 * t103 + (-t3 * t59 + t4 * t61 + (-t10 * t59 - t61 * t8) * qJD(4)) * t29; 0, qJDD(2) + t90, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t32, t33, -t32, 0, t81 * qJD(4) + t3 * t61 + t4 * t59 + t90; 0, 0, qJDD(3), t74 + t99 + t132, -t129 + t70, t52 * qJDD(3) + t59 * t83, -0.2e1 * t107 * t91 + 0.2e1 * t59 * t93, t32, t33, 0, t71 * t59 + (-t77 - t122) * t61 + t87, t71 * t61 + (t67 + t77) * t59 + t108, -t34 * qJDD(4) + (-t86 - t122) * t61 + (t23 + (t11 + (-t44 - t120) * qJD(3)) * t59) * qJD(4) + t87, -t35 * qJDD(4) + (t11 * t61 - t22 + (-t44 * t61 + t127) * qJD(3)) * qJD(4) + (t67 + t86) * t59 + t108, (-qJD(4) * t8 + qJDD(3) * t35 + t4) * t61 + (-t10 * qJD(4) + qJDD(3) * t34 - t3) * t59 + (t22 * t61 - t23 * t59 - t106 * t24 + (t34 * t61 - t35 * t59) * qJD(4)) * qJD(3) - t70, t4 * t35 - t3 * t34 - t5 * t44 - g(3) * (t109 * t45 + t46 * t44) + (t59 * t24 + t23) * t8 + (pkin(4) * t97 - t25) * t11 + (t22 - t113) * t10 + t82 * (-t109 * t46 + t44 * t45); 0, 0, 0, 0, 0, -t59 * t112, t107 * t64, t47, t93, qJDD(4), (-t6 - t100) * t59 + t75, (-t59 * t17 + t49) * qJD(4) + (t84 - t100) * t61 + t68, 0.2e1 * t101 + (t10 - t79) * qJD(4) + (pkin(4) * t112 + t65) * t59 + t75, -t64 * t127 + (t59 * t96 + t9) * qJD(4) + t65 * t61 + t68, -pkin(4) * t47 + (-t104 + t128) * t61 * qJD(3), t128 * t10 + (t3 + t90 * t61 + (-t11 * qJD(3) + t70) * t59) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t89 - t93, t47 + t83, -t106 * t64, -t81 * qJD(3) - t132 + t5;];
tau_reg = t7;
