% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RPRR2
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPRR2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR2_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR2_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR2_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:48:15
% EndTime: 2019-12-31 16:48:16
% DurationCPUTime: 0.52s
% Computational Cost: add. (764->120), mult. (1457->159), div. (0->0), fcn. (824->12), ass. (0->92)
t55 = qJ(1) + pkin(7);
t50 = qJ(3) + t55;
t45 = cos(t50);
t114 = g(2) * t45;
t53 = qJDD(1) + qJDD(3);
t113 = t53 * pkin(3);
t60 = cos(pkin(7));
t46 = t60 * pkin(1) + pkin(2);
t33 = t46 * qJD(1);
t62 = sin(qJ(3));
t106 = t62 * t33;
t59 = sin(pkin(7));
t115 = pkin(1) * t59;
t31 = t46 * qJDD(1);
t65 = cos(qJ(3));
t96 = qJD(3) * t65;
t70 = -(qJD(1) * t96 + qJDD(1) * t62) * t115 - qJD(3) * t106 + t65 * t31;
t6 = -t113 - t70;
t116 = t6 + t114;
t92 = qJD(1) * t115;
t86 = t62 * t92;
t97 = pkin(1) * qJDD(1);
t91 = t59 * t97;
t87 = -qJD(3) * t86 + t62 * t31 + t33 * t96 + t65 * t91;
t5 = t53 * pkin(6) + t87;
t61 = sin(qJ(4));
t64 = cos(qJ(4));
t15 = t65 * t92 + t106;
t54 = qJD(1) + qJD(3);
t13 = t54 * pkin(6) + t15;
t7 = t64 * qJD(2) - t61 * t13;
t98 = t7 * qJD(4);
t2 = t61 * qJDD(2) + t64 * t5 + t98;
t49 = t64 * qJDD(2);
t105 = t64 * t13;
t8 = t61 * qJD(2) + t105;
t3 = -t8 * qJD(4) - t61 * t5 + t49;
t72 = t2 * t64 - t3 * t61 + (-t61 * t8 - t64 * t7) * qJD(4);
t44 = sin(t50);
t102 = -g(1) * t45 - g(2) * t44;
t23 = t65 * t115 + t62 * t46;
t40 = g(1) * t44;
t112 = t54 * pkin(3);
t14 = t65 * t33 - t86;
t12 = -t14 - t112;
t111 = t12 * qJD(4) * t61 + t64 * t40;
t110 = t14 * t54;
t109 = t15 * t54;
t22 = -t62 * t115 + t65 * t46;
t16 = t22 * qJD(3);
t108 = t16 * t54;
t17 = t23 * qJD(3);
t107 = t17 * t54;
t104 = t64 * t53;
t103 = t45 * pkin(3) + t44 * pkin(6);
t48 = cos(t55);
t66 = cos(qJ(1));
t101 = t66 * pkin(1) + pkin(2) * t48;
t56 = t61 ^ 2;
t57 = t64 ^ 2;
t100 = t56 - t57;
t99 = t56 + t57;
t95 = qJD(4) * t64;
t58 = qJDD(2) - g(3);
t94 = t116 * t61 + t12 * t95;
t52 = t54 ^ 2;
t93 = t61 * t52 * t64;
t89 = -t44 * pkin(3) + t45 * pkin(6);
t88 = t99 * t53;
t85 = t61 * t54 * t95;
t47 = sin(t55);
t63 = sin(qJ(1));
t84 = -t63 * pkin(1) - pkin(2) * t47;
t83 = g(1) * t63 - g(2) * t66;
t81 = t7 * t61 - t8 * t64;
t80 = -t87 - t102;
t67 = qJD(4) ^ 2;
t78 = pkin(6) * t67 - t109 - t113;
t20 = -pkin(3) - t22;
t21 = pkin(6) + t23;
t77 = t20 * t53 + t21 * t67 + t107;
t76 = -pkin(6) * qJDD(4) + (t14 - t112) * qJD(4);
t75 = -qJDD(4) * t21 + (t20 * t54 - t16) * qJD(4);
t73 = -qJD(2) * qJD(4) - t12 * t54 - t102 - t5;
t71 = t102 + t72;
t69 = t40 + t70 - t114;
t28 = qJDD(4) * t64 - t67 * t61;
t27 = qJDD(4) * t61 + t67 * t64;
t19 = t57 * t53 - 0.2e1 * t85;
t18 = t56 * t53 + 0.2e1 * t85;
t11 = -0.2e1 * t100 * t54 * qJD(4) + 0.2e1 * t61 * t104;
t1 = [0, 0, 0, 0, 0, qJDD(1), t83, g(1) * t66 + g(2) * t63, 0, 0, 0, 0, 0, 0, 0, qJDD(1), g(1) * t47 - g(2) * t48 + 0.2e1 * t60 * t97, g(1) * t48 + g(2) * t47 - 0.2e1 * t91, 0, (t83 + (t59 ^ 2 + t60 ^ 2) * t97) * pkin(1), 0, 0, 0, 0, 0, t53, t22 * t53 - t107 + t69, -t23 * t53 - t108 + t80, 0, -g(1) * t84 - g(2) * t101 - t14 * t17 + t15 * t16 + t70 * t22 + t87 * t23, t18, t11, t27, t19, t28, 0, t75 * t61 + (-t77 - t116) * t64 + t111, t75 * t64 + (t77 - t40) * t61 + t94, t99 * t108 + t21 * t88 + t71, t6 * t20 + t12 * t17 - g(1) * (t84 + t89) - g(2) * (t101 + t103) - t81 * t16 + t72 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, 0, 0, 0, 0, 0, t28, -t27, 0, -t81 * qJD(4) + t2 * t61 + t3 * t64 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, t69 + t109, t80 + t110, 0, 0, t18, t11, t27, t19, t28, 0, t76 * t61 + (-t78 - t116) * t64 + t111, t76 * t64 + (t78 - t40) * t61 + t94, pkin(6) * t88 - t99 * t110 + t71, -t6 * pkin(3) + t72 * pkin(6) - g(1) * t89 - g(2) * t103 - t12 * t15 + t81 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93, t100 * t52, t61 * t53, t93, t104, qJDD(4), -g(3) * t64 + t49 + (t8 - t105) * qJD(4) + t73 * t61, t98 + (qJD(4) * t13 - t58) * t61 + t73 * t64, 0, 0;];
tau_reg = t1;
