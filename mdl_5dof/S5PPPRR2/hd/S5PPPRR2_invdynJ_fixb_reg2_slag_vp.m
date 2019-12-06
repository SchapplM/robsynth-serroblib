% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PPPRR2
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
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PPPRR2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR2_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPPRR2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:59:43
% EndTime: 2019-12-05 14:59:44
% DurationCPUTime: 0.71s
% Computational Cost: add. (843->167), mult. (1891->255), div. (0->0), fcn. (1716->10), ass. (0->99)
t51 = sin(pkin(9));
t54 = cos(pkin(9));
t52 = sin(pkin(8));
t93 = qJDD(1) * t52;
t36 = t51 * qJDD(2) + t54 * t93;
t55 = cos(pkin(8));
t43 = -t55 * qJDD(1) + qJDD(3);
t58 = sin(qJ(4));
t60 = cos(qJ(4));
t80 = t58 * t36 - t60 * t43;
t99 = qJD(1) * t52;
t38 = t51 * qJD(2) + t54 * t99;
t44 = -t55 * qJD(1) + qJD(3);
t22 = t60 * t38 + t58 * t44;
t94 = t22 * qJD(4);
t5 = -t80 - t94;
t113 = t51 * t52;
t112 = t51 * t60;
t111 = t52 * t58;
t110 = t52 * t60;
t53 = sin(pkin(7));
t109 = t53 * t55;
t56 = cos(pkin(7));
t108 = t56 * t51;
t107 = t56 * t54;
t106 = t58 * t38;
t105 = t60 * t44;
t62 = qJD(4) ^ 2;
t104 = t62 * t60;
t57 = sin(qJ(5));
t49 = t57 ^ 2;
t59 = cos(qJ(5));
t50 = t59 ^ 2;
t103 = t49 - t50;
t102 = t49 + t50;
t61 = qJD(5) ^ 2;
t101 = t61 + t62;
t100 = qJD(4) * pkin(4);
t98 = qJD(4) * t58;
t31 = t54 * t111 + t55 * t60;
t97 = qJD(5) * t31;
t96 = qJDD(4) * pkin(4);
t21 = t105 - t106;
t11 = -t21 - t100;
t95 = t11 * qJD(4);
t92 = qJDD(5) * t57;
t91 = t57 * qJDD(4);
t90 = t59 * qJDD(4);
t89 = t60 * qJDD(4);
t88 = qJD(4) * qJD(5);
t87 = t54 * t110;
t86 = t57 * t62 * t59;
t85 = -qJD(4) * t105 - t60 * t36 - t58 * t43;
t84 = t51 * t98;
t83 = -g(1) * t53 + g(2) * t56;
t82 = t57 * t88;
t81 = t59 * t88;
t79 = t102 * qJDD(4);
t78 = t57 * t81;
t28 = t54 * t109 - t108;
t30 = t55 * t107 + t53 * t51;
t77 = g(1) * (t56 * t111 + t30 * t60) + g(2) * (t53 * t111 + t28 * t60);
t76 = -g(1) * (t55 * t108 - t53 * t54) - g(2) * (t51 * t109 + t107);
t12 = qJD(4) * pkin(6) + t22;
t37 = -t54 * qJD(2) + t51 * t99;
t10 = t59 * t12 + t57 * t37;
t9 = -t57 * t12 + t59 * t37;
t75 = t10 * t59 - t9 * t57;
t74 = -t62 * t58 + t89;
t73 = qJDD(4) * t58 + t104;
t35 = -t54 * qJDD(2) + t51 * t93;
t72 = -t35 * t54 + t83;
t32 = -t55 * t58 + t87;
t17 = t59 * t113 - t32 * t57;
t18 = t57 * t113 + t32 * t59;
t34 = t59 * t112 - t57 * t54;
t33 = -t57 * t112 - t59 * t54;
t6 = -t38 * t98 - t85;
t71 = -g(1) * (t56 * t110 - t30 * t58) - g(2) * (t53 * t110 - t28 * t58) + g(3) * t31;
t70 = g(3) * t32 + t77;
t69 = g(3) * t55 + (-g(1) * t56 - g(2) * t53) * t52;
t3 = -t96 - t5;
t68 = -t3 + t71;
t4 = qJDD(4) * pkin(6) + t6;
t67 = -t4 + t77 - t95;
t66 = -pkin(6) * qJDD(5) + (t11 + t21 - t100) * qJD(5);
t1 = t9 * qJD(5) + t57 * t35 + t59 * t4;
t25 = t59 * t35;
t2 = -t10 * qJD(5) - t57 * t4 + t25;
t65 = t1 * t59 + (-t10 * t57 - t59 * t9) * qJD(5) - t2 * t57;
t64 = -pkin(6) * t61 + t68 + t94 + t96;
t63 = t65 - t70;
t24 = t31 * qJD(4);
t23 = qJD(4) * t87 - t55 * t98;
t20 = -t34 * qJD(5) + t57 * t84;
t19 = t33 * qJD(5) - t59 * t84;
t8 = t17 * qJD(5) - t24 * t59;
t7 = -t18 * qJD(5) + t24 * t57;
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1) - g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) + (t52 ^ 2 + t55 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43 * t55 - g(3) + (t35 * t51 + t36 * t54) * t52, 0, 0, 0, 0, 0, 0, -t23 * qJD(4) - t31 * qJDD(4), t24 * qJD(4) - t32 * qJDD(4), 0, t35 * t113 - t21 * t23 - t22 * t24 - t5 * t31 + t6 * t32 - g(3), 0, 0, 0, 0, 0, 0, -t31 * t90 + t7 * qJD(5) + t17 * qJDD(5) + (-t23 * t59 + t57 * t97) * qJD(4), t31 * t91 - t8 * qJD(5) - t18 * qJDD(5) + (t23 * t57 + t59 * t97) * qJD(4), (-t17 * t57 + t18 * t59) * qJDD(4) + (-t57 * t7 + t59 * t8 + (-t17 * t59 - t18 * t57) * qJD(5)) * qJD(4), t1 * t18 + t10 * t8 + t11 * t23 + t2 * t17 + t3 * t31 + t9 * t7 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) + t83, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36 * t51 + t72, 0, 0, 0, 0, 0, 0, -t73 * t51, -t74 * t51, 0, (-t5 * t58 + t6 * t60 + (-t21 * t60 - t22 * t58) * qJD(4)) * t51 + t72, 0, 0, 0, 0, 0, 0, t20 * qJD(5) + t33 * qJDD(5) + (-t59 * t104 + (t82 - t90) * t58) * t51, -t19 * qJD(5) - t34 * qJDD(5) + (t73 * t57 + t58 * t81) * t51, (-t33 * t57 + t34 * t59) * qJDD(4) + (t19 * t59 - t20 * t57 + (-t33 * t59 - t34 * t57) * qJD(5)) * qJD(4), t1 * t34 + t10 * t19 + t2 * t33 + t9 * t20 + (t3 * t58 + t60 * t95) * t51 + t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43 + t69, 0, 0, 0, 0, 0, 0, t74, -t73, 0, t5 * t60 + t6 * t58 + (-t21 * t58 + t22 * t60) * qJD(4) + t69, 0, 0, 0, 0, 0, 0, (-0.2e1 * t82 + t90) * t60 + (-t101 * t59 - t92) * t58, (-qJDD(5) * t58 - 0.2e1 * t60 * t88) * t59 + (t101 * t58 - t89) * t57, t102 * t104 + t58 * t79, (t75 * qJD(4) - t3) * t60 + (t65 + t95) * t58 + t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(4), t71 - t80, (t21 + t106) * qJD(4) + t70 + t85, 0, 0, t49 * qJDD(4) + 0.2e1 * t78, -0.2e1 * t103 * t88 + 0.2e1 * t57 * t90, t61 * t59 + t92, t50 * qJDD(4) - 0.2e1 * t78, qJDD(5) * t59 - t61 * t57, 0, t66 * t57 + t64 * t59, -t64 * t57 + t66 * t59, -t102 * t21 * qJD(4) + pkin(6) * t79 + t63, t68 * pkin(4) + t63 * pkin(6) - t11 * t22 - t75 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, t103 * t62, t91, t86, t90, qJDD(5), -g(3) * t17 + t67 * t57 + t76 * t59 + t25, g(3) * t18 + (-t35 - t76) * t57 + t67 * t59, 0, 0;];
tau_reg = t13;
