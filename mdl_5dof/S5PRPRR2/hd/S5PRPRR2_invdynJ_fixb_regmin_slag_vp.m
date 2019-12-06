% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRPRR2
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tau_reg [5x15]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRPRR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:45:15
% EndTime: 2019-12-05 15:45:17
% DurationCPUTime: 0.50s
% Computational Cost: add. (702->107), mult. (1299->156), div. (0->0), fcn. (1027->12), ass. (0->82)
t56 = qJ(2) + pkin(9) + qJ(4);
t52 = sin(t56);
t63 = sin(pkin(8));
t65 = cos(pkin(8));
t90 = g(1) * t65 + g(2) * t63;
t122 = t90 * t52;
t68 = sin(qJ(2));
t101 = qJD(1) * t68;
t71 = cos(qJ(2));
t49 = qJD(2) * pkin(2) + t71 * qJD(1);
t62 = sin(pkin(9));
t64 = cos(pkin(9));
t23 = -t62 * t101 + t64 * t49;
t22 = qJD(2) * pkin(3) + t23;
t24 = t64 * t101 + t62 * t49;
t67 = sin(qJ(4));
t70 = cos(qJ(4));
t11 = t67 * t22 + t70 * t24;
t55 = t71 * qJDD(1);
t98 = qJD(1) * qJD(2);
t35 = qJDD(2) * pkin(2) - t68 * t98 + t55;
t83 = t68 * qJDD(1) + t71 * t98;
t16 = t64 * t35 - t83 * t62;
t15 = qJDD(2) * pkin(3) + t16;
t17 = t62 * t35 + t83 * t64;
t53 = cos(t56);
t121 = g(3) * t53 + t11 * qJD(4) - t70 * t15 + t67 * t17;
t107 = t67 * t24;
t76 = g(3) * t52 - (qJD(4) * t22 + t17) * t70 + qJD(4) * t107 - t67 * t15 + t90 * t53;
t58 = qJDD(2) + qJDD(4);
t110 = t58 * pkin(4);
t119 = -t110 + t121;
t114 = pkin(2) * t62;
t54 = t64 * pkin(2) + pkin(3);
t86 = -t67 * t114 + t70 * t54;
t29 = -pkin(4) - t86;
t103 = t70 * t114 + t67 * t54;
t30 = pkin(7) + t103;
t72 = qJD(5) ^ 2;
t37 = t62 * t71 + t64 * t68;
t32 = t37 * qJD(1);
t36 = -t62 * t68 + t64 * t71;
t34 = t36 * qJD(1);
t59 = qJD(2) + qJD(4);
t94 = (-t103 * qJD(4) + t70 * t32 + t67 * t34) * t59;
t117 = t29 * t58 + t30 * t72 - t94;
t109 = t59 * pkin(4);
t108 = t11 * t59;
t69 = cos(qJ(5));
t106 = t69 * t58;
t105 = -t86 * qJD(4) - t67 * t32 + t70 * t34;
t66 = sin(qJ(5));
t60 = t66 ^ 2;
t102 = -t69 ^ 2 + t60;
t100 = t69 * qJD(5);
t99 = qJDD(1) - g(3);
t10 = t70 * t22 - t107;
t8 = -t10 - t109;
t97 = t8 * t100 + t119 * t66;
t96 = t8 * qJD(5) * t66 + t122 * t69;
t19 = t67 * t36 + t70 * t37;
t31 = t37 * qJD(2);
t33 = t36 * qJD(2);
t87 = t70 * t36 - t67 * t37;
t89 = t87 * t58 - (t19 * qJD(4) + t70 * t31 + t67 * t33) * t59;
t85 = g(1) * t63 - g(2) * t65 - qJDD(3);
t82 = pkin(7) * t72 - t108 - t110;
t81 = t19 * t72 - t89;
t80 = -pkin(7) * qJDD(5) + (t10 - t109) * qJD(5);
t4 = t87 * qJD(4) - t67 * t31 + t70 * t33;
t79 = -qJDD(5) * t19 + (-t59 * t87 - t4) * qJD(5);
t78 = -g(3) * t71 + t90 * t68;
t77 = -qJDD(5) * t30 + (t29 * t59 + t105) * qJD(5);
t75 = -t58 * pkin(7) - t8 * t59 + t76;
t74 = t122 - t121;
t73 = qJD(2) ^ 2;
t57 = t59 ^ 2;
t46 = qJDD(5) * t69 - t72 * t66;
t45 = qJDD(5) * t66 + t72 * t69;
t27 = 0.2e1 * t66 * t59 * t100 + t60 * t58;
t20 = -0.2e1 * t102 * t59 * qJD(5) + 0.2e1 * t66 * t106;
t1 = [t99, 0, t71 * qJDD(2) - t73 * t68, -qJDD(2) * t68 - t73 * t71, t16 * t36 + t17 * t37 - t23 * t31 + t24 * t33 - g(3), 0, t89, -t19 * t58 - t4 * t59, 0, 0, 0, 0, 0, t79 * t66 - t81 * t69, t81 * t66 + t79 * t69; 0, qJDD(2), t55 + t78, -t99 * t68 + t90 * t71, t23 * t32 - t24 * t34 + (t16 * t64 + t17 * t62 + t78) * pkin(2), t58, t86 * t58 + t74 + t94, -t103 * t58 + t105 * t59 + t76, t27, t20, t45, t46, 0, t77 * t66 + (-t119 - t117) * t69 + t96, t77 * t69 + (-t122 + t117) * t66 + t97; 0, 0, 0, 0, -t85, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t45; 0, 0, 0, 0, 0, t58, t74 + t108, t10 * t59 + t76, t27, t20, t45, t46, 0, t80 * t66 + (-t82 - t119) * t69 + t96, t80 * t69 + (-t122 + t82) * t66 + t97; 0, 0, 0, 0, 0, 0, 0, 0, -t66 * t57 * t69, t102 * t57, t66 * t58, t106, qJDD(5), t75 * t66 - t85 * t69, t85 * t66 + t75 * t69;];
tau_reg = t1;
