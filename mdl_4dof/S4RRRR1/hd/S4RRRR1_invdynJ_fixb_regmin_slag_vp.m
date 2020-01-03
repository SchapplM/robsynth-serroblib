% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RRRR1
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% tau_reg [4x16]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRRR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:22:15
% EndTime: 2019-12-31 17:22:17
% DurationCPUTime: 0.47s
% Computational Cost: add. (595->118), mult. (908->158), div. (0->0), fcn. (505->12), ass. (0->90)
t49 = qJ(1) + qJ(2);
t44 = qJ(3) + t49;
t33 = sin(t44);
t34 = cos(t44);
t51 = sin(qJ(3));
t52 = sin(qJ(2));
t91 = qJD(3) * t52;
t81 = qJD(1) * t91;
t77 = pkin(1) * t81;
t116 = g(1) * t34 + g(2) * t33 + t51 * t77;
t56 = cos(qJ(2));
t105 = t56 * pkin(1);
t38 = qJDD(1) * t105;
t45 = qJDD(1) + qJDD(2);
t95 = pkin(1) * qJD(1);
t86 = t52 * t95;
t14 = t45 * pkin(2) - qJD(2) * t86 + t38;
t46 = qJD(1) + qJD(2);
t19 = t46 * pkin(2) + t56 * t95;
t55 = cos(qJ(3));
t94 = qJD(2) * t56;
t82 = qJD(1) * t94;
t88 = qJDD(1) * t52;
t64 = (t82 + t88) * pkin(1);
t59 = -(qJD(3) * t19 + t64) * t55 - t51 * t14 + t116;
t90 = qJD(3) * t55;
t85 = t52 * t90;
t119 = ((t51 * t94 + t85) * qJD(1) + t51 * t88) * pkin(1);
t40 = qJDD(3) + t45;
t108 = t40 * pkin(3);
t109 = g(2) * t34;
t12 = t55 * t14;
t92 = qJD(3) * t51;
t74 = t19 * t92 - t12;
t118 = -t108 + t74 + t119 + t109;
t27 = g(1) * t33;
t114 = t27 - t109;
t101 = t52 * t55;
t72 = t51 * t56 + t101;
t17 = t72 * t95;
t35 = t51 * pkin(2) + pkin(7);
t36 = -t55 * pkin(2) - pkin(3);
t41 = qJD(3) + t46;
t58 = qJD(4) ^ 2;
t112 = t35 * t58 + t36 * t40 + (pkin(2) * t92 - t17) * t41;
t110 = pkin(2) * t40;
t107 = t41 * pkin(3);
t37 = pkin(2) + t105;
t106 = (t37 * t92 + (t72 * qJD(2) + t85) * pkin(1)) * t41;
t50 = sin(qJ(4));
t54 = cos(qJ(4));
t10 = t55 * t19 - t51 * t86;
t9 = -t10 - t107;
t104 = t9 * qJD(4) * t50 + t54 * t27;
t103 = (t51 * t19 + t55 * t86) * t41;
t102 = t51 * t52;
t100 = t54 * t40;
t99 = t55 * t37;
t98 = pkin(1) * t101 + t51 * t37;
t42 = sin(t49);
t43 = cos(t49);
t97 = g(1) * t43 + g(2) * t42;
t47 = t50 ^ 2;
t96 = -t54 ^ 2 + t47;
t89 = qJD(4) * t54;
t87 = t118 * t50 + t9 * t89;
t79 = qJD(1) * (-qJD(2) + t46);
t78 = qJD(2) * (-qJD(1) - t46);
t76 = g(1) * t42 - g(2) * t43 + t38;
t71 = t55 * t56 - t102;
t69 = pkin(7) * t58 - t103 - t108;
t15 = pkin(1) * t102 - pkin(3) - t99;
t16 = pkin(7) + t98;
t68 = t15 * t40 + t16 * t58 + t106;
t66 = -t40 * pkin(7) - t41 * t9 + t59;
t65 = -t74 + t114;
t63 = -pkin(7) * qJDD(4) + (t10 - t107) * qJD(4);
t4 = t37 * t90 + (t71 * qJD(2) - t51 * t91) * pkin(1);
t62 = -qJDD(4) * t16 + (t15 * t41 - t4) * qJD(4);
t18 = t71 * t95;
t61 = -qJDD(4) * t35 + (-pkin(2) * t90 + t36 * t41 + t18) * qJD(4);
t60 = (-pkin(2) * t41 - t19) * qJD(3) - t64;
t57 = cos(qJ(1));
t53 = sin(qJ(1));
t39 = t41 ^ 2;
t21 = qJDD(4) * t54 - t58 * t50;
t20 = qJDD(4) * t50 + t58 * t54;
t13 = 0.2e1 * t50 * t41 * t89 + t47 * t40;
t8 = -0.2e1 * t96 * t41 * qJD(4) + 0.2e1 * t50 * t100;
t1 = [qJDD(1), g(1) * t53 - g(2) * t57, g(1) * t57 + g(2) * t53, t45, (t45 * t56 + t52 * t78) * pkin(1) + t76, ((-qJDD(1) - t45) * t52 + t56 * t78) * pkin(1) + t97, t40, t40 * t99 - t106 + (-t55 * t81 + (-t82 + (-qJDD(1) - t40) * t52) * t51) * pkin(1) + t65, -t4 * t41 - t98 * t40 + t59, t13, t8, t20, t21, 0, t62 * t50 + (-t68 - t118) * t54 + t104, t62 * t54 + (t68 - t27) * t50 + t87; 0, 0, 0, t45, t52 * pkin(1) * t79 + t76, (t56 * t79 - t88) * pkin(1) + t97, t40, t17 * t41 + t12 + (-t77 + t110) * t55 + t60 * t51 + t114, t18 * t41 + (-t14 - t110) * t51 + t60 * t55 + t116, t13, t8, t20, t21, 0, t61 * t50 + (-t112 - t118) * t54 + t104, t61 * t54 + (t112 - t27) * t50 + t87; 0, 0, 0, 0, 0, 0, t40, t103 + t65 - t119, t10 * t41 + t59, t13, t8, t20, t21, 0, t63 * t50 + (-t69 - t118) * t54 + t104, t63 * t54 + (t69 - t27) * t50 + t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50 * t39 * t54, t96 * t39, t50 * t40, t100, qJDD(4), -g(3) * t54 + t66 * t50, g(3) * t50 + t66 * t54;];
tau_reg = t1;
