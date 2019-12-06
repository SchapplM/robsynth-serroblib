% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRRRR2
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
% 
% Output:
% tau_reg [5x17]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRRR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:04:53
% EndTime: 2019-12-05 17:04:55
% DurationCPUTime: 0.46s
% Computational Cost: add. (570->109), mult. (938->145), div. (0->0), fcn. (531->12), ass. (0->84)
t46 = qJD(2) + qJD(3);
t56 = cos(qJ(3));
t95 = pkin(2) * qJD(2);
t20 = pkin(3) * t46 + t56 * t95;
t51 = sin(qJ(4));
t55 = cos(qJ(4));
t52 = sin(qJ(3));
t85 = t52 * t95;
t10 = -t20 * t55 + t51 * t85;
t41 = qJD(4) + t46;
t49 = qJ(2) + qJ(3);
t44 = qJ(4) + t49;
t34 = sin(t44);
t35 = cos(t44);
t78 = qJD(4) * t85;
t114 = g(1) * t35 + g(2) * t34 + t51 * t78;
t108 = pkin(2) * t56;
t38 = qJDD(2) * t108;
t45 = qJDD(2) + qJDD(3);
t14 = pkin(3) * t45 - qJD(3) * t85 + t38;
t87 = qJDD(2) * t52;
t94 = qJD(3) * t56;
t64 = (qJD(2) * t94 + t87) * pkin(2);
t60 = -(qJD(4) * t20 + t64) * t55 - t51 * t14 + t114;
t117 = -t10 * t41 + t60;
t106 = g(2) * t35;
t12 = t55 * t14;
t91 = qJD(4) * t55;
t84 = t52 * t91;
t92 = qJD(4) * t51;
t3 = ((t51 * t94 + t84) * qJD(2) + t51 * t87) * pkin(2) + t20 * t92 - t12;
t116 = t3 + t106;
t28 = g(1) * t34;
t112 = t28 - t106;
t100 = t52 * t55;
t72 = t51 * t56 + t100;
t102 = t72 * t95 * t41;
t36 = pkin(3) * t51 + pkin(6);
t40 = qJDD(4) + t45;
t58 = qJD(5) ^ 2;
t110 = (-t40 * t55 + t41 * t92) * pkin(3) + t36 * t58 - t102;
t107 = pkin(3) * t40;
t50 = sin(qJ(5));
t54 = cos(qJ(5));
t90 = qJD(5) * t10;
t105 = t28 * t54 + t50 * t90;
t103 = (t20 * t51 + t55 * t85) * t41;
t101 = t51 * t52;
t99 = t54 * t40;
t37 = pkin(3) + t108;
t98 = pkin(2) * t100 + t37 * t51;
t42 = sin(t49);
t43 = cos(t49);
t97 = g(1) * t43 + g(2) * t42;
t47 = t50 ^ 2;
t96 = -t54 ^ 2 + t47;
t89 = qJD(5) * t54;
t88 = qJDD(1) - g(3);
t86 = t10 * t89 + t116 * t50;
t80 = qJD(2) * (-qJD(3) + t46);
t79 = qJD(3) * (-qJD(2) - t46);
t77 = g(1) * t42 - g(2) * t43 + t38;
t75 = pkin(6) * t58 - t103;
t18 = pkin(2) * t101 - t37 * t55;
t74 = -t18 * t40 - (t37 * t92 + (qJD(3) * t72 + t84) * pkin(2)) * t41;
t71 = t55 * t56 - t101;
t70 = -pkin(6) * qJDD(5) - t90;
t15 = pkin(6) + t98;
t67 = t15 * t58 - t74;
t65 = -t40 * pkin(6) + t117;
t4 = t37 * t91 + (qJD(3) * t71 - t52 * t92) * pkin(2);
t63 = -qJDD(5) * t15 + (t18 * t41 - t4) * qJD(5);
t17 = t71 * t95;
t62 = -qJDD(5) * t36 + (t17 + (-qJD(4) - t41) * t55 * pkin(3)) * qJD(5);
t61 = (-pkin(3) * t41 - t20) * qJD(4) - t64;
t59 = t112 - t3;
t57 = cos(qJ(2));
t53 = sin(qJ(2));
t39 = t41 ^ 2;
t22 = qJDD(5) * t54 - t50 * t58;
t21 = qJDD(5) * t50 + t54 * t58;
t13 = 0.2e1 * t41 * t50 * t89 + t40 * t47;
t6 = -0.2e1 * qJD(5) * t41 * t96 + 0.2e1 * t50 * t99;
t1 = [t88, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t21; 0, qJDD(2), g(1) * t53 - g(2) * t57, g(1) * t57 + g(2) * t53, t45, (t45 * t56 + t52 * t79) * pkin(2) + t77, ((-qJDD(2) - t45) * t52 + t56 * t79) * pkin(2) + t97, t40, t59 + t74, -t4 * t41 - t40 * t98 + t60, t13, t6, t21, t22, 0, t63 * t50 + (-t67 - t116) * t54 + t105, t63 * t54 + (t67 - t28) * t50 + t86; 0, 0, 0, 0, t45, pkin(2) * t52 * t80 + t77, (t56 * t80 - t87) * pkin(2) + t97, t40, t102 + t12 + (-t78 + t107) * t55 + t61 * t51 + t112, t17 * t41 + (-t14 - t107) * t51 + t61 * t55 + t114, t13, t6, t21, t22, 0, t62 * t50 + (-t110 - t116) * t54 + t105, t62 * t54 + (t110 - t28) * t50 + t86; 0, 0, 0, 0, 0, 0, 0, t40, t59 + t103, t117, t13, t6, t21, t22, 0, t70 * t50 + (-t75 - t116) * t54 + t105, t70 * t54 + (t75 - t28) * t50 + t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50 * t39 * t54, t96 * t39, t50 * t40, t99, qJDD(5), t50 * t65 + t54 * t88, -t50 * t88 + t54 * t65;];
tau_reg = t1;
