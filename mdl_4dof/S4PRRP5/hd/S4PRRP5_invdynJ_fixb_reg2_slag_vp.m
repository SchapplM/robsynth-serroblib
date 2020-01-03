% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4PRRP5
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
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PRRP5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP5_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:29:25
% EndTime: 2019-12-31 16:29:26
% DurationCPUTime: 0.65s
% Computational Cost: add. (435->151), mult. (1037->204), div. (0->0), fcn. (617->6), ass. (0->99)
t49 = sin(qJ(2));
t108 = g(3) * t49;
t51 = cos(qJ(2));
t45 = sin(pkin(6));
t46 = cos(pkin(6));
t68 = g(1) * t46 + g(2) * t45;
t60 = t68 * t51;
t55 = -t108 - t60;
t61 = t68 * t49;
t89 = qJDD(1) - g(3);
t118 = t51 * t89 + t61;
t32 = qJD(2) * pkin(5) + t49 * qJD(1);
t90 = t51 * qJD(1);
t95 = qJD(2) * pkin(2);
t33 = -t90 - t95;
t48 = sin(qJ(3));
t43 = t48 ^ 2;
t50 = cos(qJ(3));
t44 = t50 ^ 2;
t97 = t43 + t44;
t81 = t97 * t51;
t117 = t32 * t81 + t33 * t49;
t71 = qJ(4) * qJD(2) + t32;
t9 = t71 * t48;
t10 = t71 * t50;
t42 = t50 * qJDD(2);
t83 = qJD(2) * qJD(3);
t79 = t48 * t83;
t114 = -0.2e1 * t79 + t42;
t94 = qJD(3) * pkin(3);
t8 = -t9 + t94;
t67 = t10 * t50 - t48 * t8;
t113 = qJD(2) * t67;
t112 = t9 + t8;
t111 = pkin(3) * t43;
t107 = g(3) * t51;
t106 = t50 * pkin(3);
t104 = t48 * t51;
t103 = t49 * t50;
t102 = t50 * t51;
t53 = qJD(2) ^ 2;
t101 = t50 * t53;
t100 = qJ(4) + pkin(5);
t80 = qJD(3) * t90;
t99 = g(3) * t104 + t50 * t80;
t98 = t43 - t44;
t52 = qJD(3) ^ 2;
t96 = t52 + t53;
t93 = qJDD(2) * pkin(2);
t92 = qJDD(3) * pkin(3);
t40 = pkin(2) + t106;
t17 = -qJD(2) * t40 + qJD(4) - t90;
t91 = t17 * qJD(2);
t88 = qJDD(2) * t40;
t87 = qJDD(3) * t48;
t41 = t48 * qJDD(2);
t86 = t51 * qJDD(1);
t85 = t51 * qJDD(2);
t84 = qJD(1) * qJD(2);
t19 = qJDD(2) * pkin(5) + t49 * qJDD(1) + t51 * t84;
t82 = t97 * t19;
t39 = t49 * t84;
t78 = t50 * t83;
t77 = qJD(1) * t97;
t76 = qJD(3) * t100;
t75 = t68 * t103 + t50 * t39 + t48 * t80;
t74 = -qJD(2) * t33 - t19;
t56 = pkin(3) * t79 + qJDD(4) + t39 - t88;
t6 = t56 - t86;
t73 = t6 - t88;
t72 = t97 * qJDD(2);
t70 = -qJ(4) * qJDD(2) - t19;
t69 = t48 * t78;
t65 = -g(1) * (-t104 * t46 + t45 * t50) - g(2) * (-t104 * t45 - t46 * t50) + t48 * t108;
t64 = -g(1) * (-t102 * t46 - t45 * t48) - g(2) * (-t102 * t45 + t46 * t48) + g(3) * t103;
t63 = qJD(3) * t71;
t18 = t39 - t86 - t93;
t62 = pkin(5) * t52 + t18 - t93;
t59 = qJD(2) * qJD(4) - t70;
t58 = (-qJD(4) - t17) * qJD(2) + t70;
t57 = -pkin(5) * qJDD(3) + (t33 - t95) * qJD(3);
t54 = (-t68 - t84) * t49;
t34 = t48 * t101;
t28 = t100 * t50;
t27 = t100 * t48;
t24 = t98 * t53;
t23 = qJDD(3) * t50 - t52 * t48;
t22 = t52 * t50 + t87;
t21 = t44 * qJDD(2) - 0.2e1 * t69;
t20 = t43 * qJDD(2) + 0.2e1 * t69;
t12 = -t48 * qJD(4) - t50 * t76;
t11 = t50 * qJD(4) - t48 * t76;
t7 = 0.2e1 * t42 * t48 - 0.2e1 * t83 * t98;
t5 = t49 * t72 + t53 * t81;
t4 = t114 * t51 + (-t50 * t96 - t87) * t49;
t3 = (-qJDD(3) * t49 - 0.2e1 * t51 * t83) * t50 + (t49 * t96 - t85) * t48;
t2 = -t48 * t63 + t50 * t59;
t1 = -t48 * t59 - t50 * t63 + t92;
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t89, 0, 0, 0, 0, 0, 0, -t53 * t49 + t85, -qJDD(2) * t49 - t53 * t51, 0, -g(3) + (t49 ^ 2 + t51 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, t4, t3, t5, t117 * qJD(2) - t18 * t51 + t49 * t82 - g(3), 0, 0, 0, 0, 0, 0, t4, t3, t5, -g(3) + (-t6 + t113) * t51 + (t91 - t1 * t48 + t2 * t50 + (-t10 * t48 - t50 * t8) * qJD(3)) * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t118, -t49 * t89 + t60, 0, 0, t20, t7, t22, t21, t23, 0, t57 * t48 + (-t62 - t107) * t50 + t75, t57 * t50 + (t54 + t62) * t48 + t99, -t108 + t82 + pkin(5) * t72 + (-qJD(2) * t77 - t68) * t51, (-t18 + t61 - t107) * pkin(2) + (t82 + t55) * pkin(5) - t117 * qJD(1), t20, t7, t22, t21, t23, 0, -t27 * qJDD(3) + (-t73 - t107) * t50 + (t12 + (t17 + (-t40 - t106) * qJD(2)) * t48) * qJD(3) + t75, -t28 * qJDD(3) + (t17 * t50 - t11 + (-t40 * t50 + t111) * qJD(2)) * qJD(3) + (t54 + t73) * t48 + t99, (-qJD(3) * t8 + qJDD(2) * t28 + t2) * t50 + (-t10 * qJD(3) + qJDD(2) * t27 - t1) * t48 + (t11 * t50 - t12 * t48 + (t27 * t50 - t28 * t48) * qJD(3) - t51 * t77) * qJD(2) + t55, t2 * t28 + t10 * t11 - t1 * t27 + t8 * t12 - t6 * t40 + t17 * t48 * t94 - g(3) * (t100 * t49 + t51 * t40) + (-t17 * t49 - t51 * t67) * qJD(1) + t68 * (-t100 * t51 + t40 * t49); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, t24, t41, t34, t42, qJDD(3), t48 * t74 + t65, t50 * t74 + t64, 0, 0, -t34, t24, t41, t34, t42, qJDD(3), 0.2e1 * t92 + (pkin(3) * t101 + t58) * t48 + t65, -t53 * t111 + t58 * t50 + t64, -pkin(3) * t41 + (-t94 + t112) * t50 * qJD(2), t112 * t10 + (t1 + (-g(1) * t45 + g(2) * t46) * t50 + (-t91 - t55) * t48) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t114, t41 + 0.2e1 * t78, -t97 * t53, -t113 + t56 - t118;];
tau_reg = t13;
