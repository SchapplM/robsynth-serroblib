% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRRR4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% tau_reg [5x17]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRR4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:52:24
% EndTime: 2020-01-03 11:52:26
% DurationCPUTime: 0.58s
% Computational Cost: add. (1051->121), mult. (1860->156), div. (0->0), fcn. (1095->14), ass. (0->91)
t60 = cos(pkin(9));
t47 = pkin(1) * t60 + pkin(2);
t39 = t47 * qJD(1);
t63 = sin(qJ(3));
t67 = cos(qJ(3));
t59 = sin(pkin(9));
t120 = pkin(1) * t59;
t98 = qJD(1) * t120;
t23 = t63 * t39 + t67 * t98;
t66 = cos(qJ(4));
t106 = t66 * t23;
t91 = t63 * t98;
t22 = t67 * t39 - t91;
t56 = qJD(1) + qJD(3);
t20 = pkin(3) * t56 + t22;
t62 = sin(qJ(4));
t11 = t62 * t20 + t106;
t55 = qJDD(1) + qJDD(3);
t103 = qJD(3) * t39;
t37 = t47 * qJDD(1);
t86 = -t63 * t103 + t67 * t37;
t96 = qJD(1) * qJD(3) * t67;
t73 = (-qJDD(1) * t63 - t96) * t120 + t86;
t13 = pkin(3) * t55 + t73;
t16 = (qJDD(1) * t120 + t103) * t67 - qJD(3) * t91 + t37 * t63;
t54 = qJ(1) + pkin(9) + qJ(3);
t48 = qJ(4) + t54;
t43 = sin(t48);
t44 = cos(t48);
t78 = -g(2) * t44 - g(3) * t43 + t66 * t13 - t62 * t16;
t71 = -t11 * qJD(4) + t78;
t52 = qJDD(4) + t55;
t118 = pkin(4) * t52;
t83 = t118 + t71;
t102 = qJD(4) * t62;
t90 = g(2) * t43 - g(3) * t44 + t23 * t102;
t74 = -(qJD(4) * t20 + t16) * t66 - t62 * t13 + t90;
t40 = t67 * t47;
t125 = -t63 * t120 + t40;
t29 = pkin(3) + t125;
t30 = t67 * t120 + t47 * t63;
t105 = t62 * t29 + t66 * t30;
t45 = sin(t54);
t46 = cos(t54);
t123 = -g(2) * t46 - g(3) * t45;
t14 = t22 * t62 + t106;
t49 = pkin(3) * t62 + pkin(8);
t119 = pkin(3) * t66;
t50 = -pkin(4) - t119;
t53 = qJD(4) + t56;
t69 = qJD(5) ^ 2;
t122 = t49 * t69 + t50 * t52 + (pkin(3) * t102 - t14) * t53;
t117 = pkin(4) * t53;
t27 = t125 * qJD(3);
t28 = t30 * qJD(3);
t111 = (t105 * qJD(4) + t27 * t62 + t28 * t66) * t53;
t110 = t11 * t53;
t109 = t23 * t62;
t65 = cos(qJ(5));
t107 = t65 * t52;
t61 = sin(qJ(5));
t57 = t61 ^ 2;
t104 = -t65 ^ 2 + t57;
t101 = t65 * qJD(5);
t100 = qJDD(2) - g(1);
t97 = -pkin(3) * t53 - t20;
t10 = t20 * t66 - t109;
t8 = -t10 - t117;
t95 = t8 * t101 - t83 * t61;
t64 = sin(qJ(1));
t68 = cos(qJ(1));
t87 = -g(2) * t68 - g(3) * t64;
t85 = t29 * t66 - t30 * t62;
t81 = pkin(8) * t69 - t110 - t118;
t17 = -pkin(4) - t85;
t18 = pkin(8) + t105;
t80 = t17 * t52 + t18 * t69 + t111;
t79 = -pkin(8) * t52 - t8 * t53 + t74;
t77 = -pkin(8) * qJDD(5) + (t10 - t117) * qJD(5);
t4 = t85 * qJD(4) + t27 * t66 - t28 * t62;
t76 = -qJDD(5) * t18 + (t17 * t53 - t4) * qJD(5);
t15 = t22 * t66 - t109;
t75 = -qJDD(5) * t49 + (-qJD(4) * t119 + t50 * t53 + t15) * qJD(5);
t72 = g(2) * t45 - g(3) * t46 - t16;
t51 = t53 ^ 2;
t36 = qJDD(5) * t65 - t61 * t69;
t35 = qJDD(5) * t61 + t65 * t69;
t24 = 0.2e1 * t53 * t61 * t101 + t52 * t57;
t19 = -0.2e1 * t104 * t53 * qJD(5) + 0.2e1 * t61 * t107;
t6 = t8 * qJD(5) * t61;
t1 = [qJDD(1), t87, g(2) * t64 - g(3) * t68, (t87 + (t59 ^ 2 + t60 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t55, -t28 * t56 + t40 * t55 + (-t96 + (-qJDD(1) - t55) * t63) * t120 + t86 + t123, -t27 * t56 - t30 * t55 + t72, t52, t85 * t52 - t111 + t71, -t105 * t52 - t4 * t53 + t74, t24, t19, t35, t36, 0, t6 + t76 * t61 + (-t80 + t83) * t65, t80 * t61 + t76 * t65 + t95; 0, 0, 0, t100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t35; 0, 0, 0, 0, t55, t23 * t56 + t123 + t73, t22 * t56 + t72, t52, t52 * t119 + t14 * t53 + (t97 * t62 - t106) * qJD(4) + t78, t15 * t53 + (-pkin(3) * t52 - t13) * t62 + (t97 * qJD(4) - t16) * t66 + t90, t24, t19, t35, t36, 0, t6 + t75 * t61 + (-t122 + t83) * t65, t122 * t61 + t75 * t65 + t95; 0, 0, 0, 0, 0, 0, 0, t52, t71 + t110, t10 * t53 + t74, t24, t19, t35, t36, 0, t6 + t77 * t61 + (-t81 + t83) * t65, t81 * t61 + t77 * t65 + t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61 * t51 * t65, t104 * t51, t61 * t52, t107, qJDD(5), t100 * t65 + t79 * t61, -t100 * t61 + t79 * t65;];
tau_reg = t1;
