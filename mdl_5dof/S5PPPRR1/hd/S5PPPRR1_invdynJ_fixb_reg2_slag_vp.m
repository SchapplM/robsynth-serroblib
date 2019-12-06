% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PPPRR1
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
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PPPRR1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR1_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPPRR1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:58:12
% EndTime: 2019-12-05 14:58:13
% DurationCPUTime: 0.62s
% Computational Cost: add. (820->132), mult. (1835->202), div. (0->0), fcn. (1661->12), ass. (0->92)
t57 = cos(pkin(8));
t43 = -t57 * qJDD(1) + qJDD(3);
t54 = sin(pkin(8));
t55 = sin(pkin(7));
t58 = cos(pkin(7));
t75 = g(3) * t57 - (g(1) * t58 + g(2) * t55) * t54;
t17 = t43 + t75;
t53 = sin(pkin(9));
t56 = cos(pkin(9));
t90 = qJDD(1) * t54;
t32 = qJDD(2) * t56 - t53 * t90;
t33 = qJDD(2) * t53 + t56 * t90;
t60 = sin(qJ(4));
t62 = cos(qJ(4));
t82 = -t62 * t32 + t60 * t33;
t98 = qJD(1) * t54;
t34 = qJD(2) * t56 - t53 * t98;
t35 = qJD(2) * t53 + t56 * t98;
t14 = t34 * t60 + t35 * t62;
t92 = t14 * qJD(4);
t6 = -t82 - t92;
t37 = t53 * t62 + t56 * t60;
t31 = t37 * qJD(4);
t107 = g(3) * t54;
t106 = t55 * t57;
t50 = pkin(9) + qJ(4);
t46 = sin(t50);
t105 = t58 * t46;
t47 = cos(t50);
t104 = t58 * t47;
t59 = sin(qJ(5));
t61 = cos(qJ(5));
t103 = t59 * t61;
t102 = t60 * t35;
t51 = t59 ^ 2;
t52 = t61 ^ 2;
t101 = t51 - t52;
t100 = t51 + t52;
t99 = qJD(4) * pkin(4);
t36 = t53 * t60 - t56 * t62;
t97 = qJD(4) * t36;
t96 = qJD(4) * t60;
t95 = qJD(4) * t62;
t24 = t37 * t54;
t94 = qJD(5) * t24;
t93 = qJDD(4) * pkin(4);
t91 = t97 * qJD(4);
t89 = t59 * qJDD(4);
t88 = t61 * qJDD(4);
t87 = qJD(4) * qJD(5);
t64 = qJD(4) ^ 2;
t86 = t64 * t103;
t85 = -t60 * t32 - t62 * t33 - t34 * t95;
t84 = -g(1) * t55 + g(2) * t58;
t83 = -t43 * t57 - g(3);
t81 = t100 * qJDD(4);
t80 = t87 * t103;
t12 = qJD(4) * pkin(6) + t14;
t44 = -qJD(1) * t57 + qJD(3);
t10 = t12 * t61 + t44 * t59;
t9 = -t12 * t59 + t44 * t61;
t78 = t10 * t61 - t59 * t9;
t25 = t36 * t54;
t16 = -t25 * t61 - t57 * t59;
t15 = t25 * t59 - t57 * t61;
t13 = t34 * t62 - t102;
t77 = -qJD(4) * t31 - qJDD(4) * t36;
t5 = -t35 * t96 - t85;
t74 = -g(1) * (-t105 * t57 + t47 * t55) - g(2) * (-t106 * t46 - t104) + t46 * t107;
t73 = g(1) * (t104 * t57 + t46 * t55) + g(2) * (t106 * t47 - t105) + t47 * t107;
t63 = qJD(5) ^ 2;
t72 = t37 * t63 - t77;
t4 = -t93 - t6;
t71 = -t4 + t74;
t70 = 0.2e1 * qJD(5) * t97 - qJDD(5) * t37;
t11 = -t13 - t99;
t69 = -pkin(6) * qJDD(5) + (t11 + t13 - t99) * qJD(5);
t3 = qJDD(4) * pkin(6) + t5;
t68 = -qJD(4) * t11 - t3 + t73;
t1 = t9 * qJD(5) + t61 * t3 + t59 * t43;
t38 = t61 * t43;
t2 = -t10 * qJD(5) - t59 * t3 + t38;
t67 = t1 * t61 - t2 * t59 + (-t10 * t59 - t61 * t9) * qJD(5);
t66 = -pkin(6) * t63 + t71 + t92 + t93;
t65 = t67 - t73;
t41 = qJDD(5) * t61 - t59 * t63;
t40 = qJDD(5) * t59 + t61 * t63;
t19 = (-t53 * t96 + t56 * t95) * t54;
t18 = t54 * t31;
t8 = -qJD(5) * t16 + t59 * t18;
t7 = qJD(5) * t15 - t61 * t18;
t20 = [0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1) - g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) + (t54 ^ 2 + t57 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t32 * t53 + t33 * t56) * t54 + t83, 0, 0, 0, 0, 0, 0, -qJD(4) * t19 - qJDD(4) * t24, qJD(4) * t18 + qJDD(4) * t25, 0, -t13 * t19 - t14 * t18 - t24 * t6 - t25 * t5 + t83, 0, 0, 0, 0, 0, 0, -t24 * t88 + t8 * qJD(5) + t15 * qJDD(5) + (-t19 * t61 + t59 * t94) * qJD(4), t24 * t89 - t7 * qJD(5) - t16 * qJDD(5) + (t19 * t59 + t61 * t94) * qJD(4), (-t15 * t59 + t16 * t61) * qJDD(4) + (-t59 * t8 + t61 * t7 + (-t15 * t61 - t16 * t59) * qJD(5)) * qJD(4), t1 * t16 + t10 * t7 + t11 * t19 + t15 * t2 + t24 * t4 + t8 * t9 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) + t84, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32 * t56 + t33 * t53 + t84, 0, 0, 0, 0, 0, 0, t77, -qJDD(4) * t37 + t91, 0, -t13 * t31 - t14 * t97 - t36 * t6 + t37 * t5 + t84, 0, 0, 0, 0, 0, 0, t59 * t70 - t61 * t72, t59 * t72 + t61 * t70, -t100 * t91 + t37 * t81, t11 * t31 + t4 * t36 + t37 * t67 - t78 * t97 + t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, 0, 0, 0, 0, 0, t41, -t40, 0, qJD(5) * t78 + t1 * t59 + t2 * t61 + t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(4), t74 - t82, (t13 + t102) * qJD(4) + t73 + t85, 0, 0, qJDD(4) * t51 + 0.2e1 * t80, -0.2e1 * t101 * t87 + 0.2e1 * t59 * t88, t40, qJDD(4) * t52 - 0.2e1 * t80, t41, 0, t59 * t69 + t61 * t66, -t59 * t66 + t61 * t69, -qJD(4) * t100 * t13 + pkin(6) * t81 + t65, pkin(4) * t71 + pkin(6) * t65 - t11 * t14 - t13 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, t101 * t64, t89, t86, t88, qJDD(5), t68 * t59 + t61 * t75 + t38, -t17 * t59 + t68 * t61, 0, 0;];
tau_reg = t20;
