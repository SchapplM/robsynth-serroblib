% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RPRP3
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
% Datum: 2019-12-31 16:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPRP3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP3_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:42:50
% EndTime: 2019-12-31 16:42:51
% DurationCPUTime: 0.53s
% Computational Cost: add. (494->149), mult. (1015->188), div. (0->0), fcn. (550->8), ass. (0->97)
t53 = sin(pkin(6));
t37 = t53 * pkin(1) + pkin(5);
t25 = t37 * qJDD(1);
t108 = qJD(2) * qJD(3) + t25;
t58 = cos(qJ(3));
t27 = t37 * qJD(1);
t75 = qJ(4) * qJD(1) + t27;
t68 = t75 * t58;
t90 = pkin(1) * qJDD(1);
t47 = t58 * qJD(2);
t56 = sin(qJ(3));
t7 = -t75 * t56 + t47;
t92 = qJD(3) * pkin(3);
t5 = t7 + t92;
t107 = t5 - t7;
t51 = t56 ^ 2;
t106 = pkin(3) * t51;
t50 = qJ(1) + pkin(6);
t41 = sin(t50);
t105 = g(1) * t41;
t42 = cos(t50);
t104 = g(1) * t42;
t103 = g(2) * t41;
t102 = g(2) * t42;
t101 = g(3) * t58;
t54 = cos(pkin(6));
t100 = t54 * pkin(1);
t57 = sin(qJ(1));
t99 = t57 * pkin(1);
t98 = t58 * pkin(3);
t97 = t41 * t58;
t96 = t42 * t56;
t61 = qJD(1) ^ 2;
t95 = t58 * t61;
t52 = t58 ^ 2;
t94 = t51 - t52;
t93 = t51 + t52;
t91 = qJ(4) + t37;
t38 = -pkin(2) - t100;
t28 = qJD(1) * t38;
t89 = qJD(1) * t58;
t88 = qJD(3) * t56;
t87 = qJDD(3) * pkin(3);
t86 = t56 * qJD(2);
t40 = pkin(2) + t98;
t21 = -t40 - t100;
t14 = qJD(1) * t21 + qJD(4);
t85 = -qJD(4) - t14;
t84 = qJDD(1) * t21;
t26 = qJDD(1) * t38;
t44 = t56 * qJDD(1);
t46 = t58 * qJDD(1);
t83 = qJ(4) * qJDD(1);
t82 = qJD(1) * qJD(3);
t81 = qJD(1) * qJD(4);
t79 = t56 * t82;
t78 = t58 * t82;
t3 = t56 * qJDD(2) + t108 * t58 - t27 * t88;
t6 = pkin(3) * t79 + qJDD(4) + t84;
t77 = t6 + t84;
t76 = qJD(3) * t91;
t74 = t56 * t78;
t73 = t103 + t104;
t72 = t102 - t105;
t59 = cos(qJ(1));
t71 = g(1) * t57 - g(2) * t59;
t8 = t86 + t68;
t70 = t5 * t56 - t58 * t8;
t45 = t58 * qJDD(2);
t69 = g(1) * t96 + t56 * t103 - t101 + t45;
t13 = t58 * t27 + t86;
t60 = qJD(3) ^ 2;
t67 = t37 * t60 + 0.2e1 * t26;
t66 = g(2) * t97 + g(3) * t56 + t58 * t104 - t3;
t65 = -t83 - t108;
t64 = 0.2e1 * qJD(3) * t28 - qJDD(3) * t37;
t12 = -t56 * t27 + t47;
t4 = -t13 * qJD(3) - t56 * t25 + t45;
t63 = t3 * t58 - t4 * t56 + (-t12 * t58 - t13 * t56) * qJD(3);
t55 = -qJ(4) - pkin(5);
t48 = t59 * pkin(1);
t35 = t56 * t95;
t32 = g(1) * t97;
t31 = g(2) * t96;
t24 = t94 * t61;
t23 = qJDD(3) * t58 - t60 * t56;
t22 = qJDD(3) * t56 + t60 * t58;
t18 = t52 * qJDD(1) - 0.2e1 * t74;
t17 = t51 * qJDD(1) + 0.2e1 * t74;
t16 = t91 * t58;
t15 = t91 * t56;
t11 = -t56 * qJD(4) - t58 * t76;
t10 = t58 * qJD(4) - t56 * t76;
t9 = 0.2e1 * t56 * t46 - 0.2e1 * t94 * t82;
t2 = t58 * t81 + (-t79 + t46) * qJ(4) + t3;
t1 = t87 + t45 - qJD(3) * t68 + (t65 - t81) * t56;
t19 = [0, 0, 0, 0, 0, qJDD(1), t71, g(1) * t59 + g(2) * t57, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t54 * t90 - t72, -0.2e1 * t53 * t90 + t73, 0, (t71 + (t53 ^ 2 + t54 ^ 2) * t90) * pkin(1), t17, t9, t22, t18, t23, 0, t32 + t64 * t56 + (-t67 - t102) * t58, t31 + t64 * t58 + (t67 - t105) * t56, t93 * t25 + t63 - t73, t26 * t38 - g(1) * (-t41 * pkin(2) + t42 * pkin(5) - t99) - g(2) * (t42 * pkin(2) + t41 * pkin(5) + t48) + t63 * t37, t17, t9, t22, t18, t23, 0, -t15 * qJDD(3) + t32 + (-t77 - t102) * t58 + (t11 + (t14 + (t21 - t98) * qJD(1)) * t56) * qJD(3), -t16 * qJDD(3) + t31 + (t77 - t105) * t56 + (t14 * t58 - t10 + (t21 * t58 + t106) * qJD(1)) * qJD(3), (-qJD(3) * t5 + qJDD(1) * t16 + t2 + (qJD(3) * t15 + t10) * qJD(1)) * t58 + (-t8 * qJD(3) + qJDD(1) * t15 - t1 + (-qJD(3) * t16 - t11) * qJD(1)) * t56 - t73, t2 * t16 + t8 * t10 - t1 * t15 + t5 * t11 + t6 * t21 + t14 * pkin(3) * t88 - g(1) * (-t41 * t40 - t42 * t55 - t99) - g(2) * (t42 * t40 - t41 * t55 + t48); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, 0, t23, -t22, 0, t3 * t56 + t4 * t58 - g(3) + (-t12 * t56 + t13 * t58) * qJD(3), 0, 0, 0, 0, 0, 0, t23, -t22, 0, -t70 * qJD(3) + t1 * t58 + t2 * t56 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, t24, t44, t35, t46, qJDD(3), (-qJD(1) * t28 - t25) * t56 + t69, t12 * qJD(3) - t28 * t89 + t66, 0, 0, -t35, t24, t44, t35, t46, qJDD(3), 0.2e1 * t87 + (t8 - t68) * qJD(3) + (pkin(3) * t95 + t85 * qJD(1) + t65) * t56 + t69, -t61 * t106 - t58 * t83 + t7 * qJD(3) + (qJ(4) * t88 + t85 * t58) * qJD(1) + t66, -pkin(3) * t44 + (-t92 + t107) * t89, t107 * t8 + (-t101 + t1 + (-qJD(1) * t14 + t73) * t56) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46 + 0.2e1 * t79, t44 + 0.2e1 * t78, -t93 * t61, t70 * qJD(1) + t6 + t72;];
tau_reg = t19;
