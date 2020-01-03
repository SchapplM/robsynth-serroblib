% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RPRP4
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
% Datum: 2019-12-31 16:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPRP4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP4_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP4_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP4_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:43:55
% EndTime: 2019-12-31 16:43:56
% DurationCPUTime: 0.46s
% Computational Cost: add. (493->121), mult. (981->157), div. (0->0), fcn. (528->8), ass. (0->81)
t106 = 2 * qJD(3);
t49 = sin(pkin(6));
t36 = t49 * pkin(1) + pkin(5);
t105 = (qJD(2) * qJD(3)) + qJDD(1) * t36;
t46 = qJ(1) + pkin(6);
t40 = sin(t46);
t41 = cos(t46);
t71 = g(1) * t41 + g(2) * t40;
t53 = cos(qJ(3));
t27 = t36 * qJD(1);
t51 = sin(qJ(3));
t96 = t51 * t27;
t10 = t53 * qJD(2) - t96;
t104 = qJD(4) - t10;
t6 = -(qJD(3) * pkin(3)) + t104;
t11 = t51 * qJD(2) + t53 * t27;
t7 = (qJD(3) * qJ(4)) + t11;
t67 = pkin(3) * t51 - qJ(4) * t53;
t12 = t67 * qJD(3) - t51 * qJD(4);
t68 = t53 * pkin(3) + t51 * qJ(4);
t66 = -pkin(2) - t68;
t50 = cos(pkin(6));
t97 = t50 * pkin(1);
t13 = t66 - t97;
t3 = qJD(1) * t12 + qJDD(1) * t13;
t93 = pkin(1) * qJDD(1);
t89 = qJDD(3) * pkin(3);
t103 = qJDD(4) - t89;
t8 = qJD(1) * t13;
t86 = qJDD(3) * t36;
t102 = t8 * t106 - t86;
t101 = g(1) * t40;
t98 = g(2) * t41;
t95 = t51 * t53;
t47 = t51 ^ 2;
t48 = t53 ^ 2;
t94 = -t47 + t48;
t37 = -pkin(2) - t97;
t28 = qJD(1) * t37;
t91 = qJD(1) * t51;
t90 = qJD(3) * t27;
t26 = qJDD(1) * t37;
t85 = t47 * qJDD(1);
t84 = t48 * qJDD(1);
t83 = t53 * qJDD(1);
t82 = qJD(1) * qJD(3);
t80 = qJDD(3) * qJ(4);
t79 = t51 * qJDD(2) + t105 * t53;
t54 = cos(qJ(1));
t78 = t54 * pkin(1) + t41 * pkin(2) + t40 * pkin(5);
t52 = sin(qJ(1));
t77 = -t52 * pkin(1) + t41 * pkin(5);
t76 = t10 + t96;
t74 = (-qJDD(2) + t90) * t53 + t105 * t51;
t72 = t82 * t95;
t70 = g(1) * t52 - g(2) * t54;
t55 = qJD(3) ^ 2;
t69 = t36 * t55 + t98;
t65 = -t71 + (t84 + t85) * t36;
t64 = -g(3) * t53 + t71 * t51 - t74;
t63 = -0.2e1 * t26 - t69;
t62 = t28 * t106 - t86;
t61 = t11 * qJD(3) + t64;
t60 = -0.2e1 * t3 - t69;
t1 = t80 + (qJD(4) - t96) * qJD(3) + t79;
t2 = t74 + t103;
t59 = t1 * t53 + t2 * t51 + (-t51 * t7 + t53 * t6) * qJD(3);
t4 = -t51 * t90 + t79;
t58 = t4 * t53 + t74 * t51 + (-t10 * t53 - t11 * t51) * qJD(3);
t56 = qJD(1) ^ 2;
t43 = t51 * qJDD(1);
t32 = t56 * t95;
t30 = t53 * t101;
t24 = t94 * t56;
t23 = qJDD(3) * t53 - t55 * t51;
t22 = qJDD(3) * t51 + t55 * t53;
t19 = t67 * qJD(1);
t15 = -0.2e1 * t72 + t84;
t14 = 0.2e1 * t72 + t85;
t9 = t51 * t83 + t94 * t82;
t5 = [0, 0, 0, 0, 0, qJDD(1), t70, g(1) * t54 + g(2) * t52, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t50 * t93 + t101 - t98, -0.2e1 * t49 * t93 + t71, 0, (t70 + (t49 ^ 2 + t50 ^ 2) * t93) * pkin(1), t14, 0.2e1 * t9, t22, t15, t23, 0, t62 * t51 + t63 * t53 + t30, t62 * t53 + (-t63 - t101) * t51, t58 + t65, t26 * t37 - g(1) * (-t40 * pkin(2) + t77) - g(2) * t78 + t58 * t36, t14, t22, -0.2e1 * t9, 0, -t23, t15, t102 * t51 + t60 * t53 + t30, t59 + t65, -t102 * t53 + (t60 + t101) * t51, t3 * t13 + t8 * t12 - g(1) * t77 - g(2) * (t68 * t41 + t78) - t66 * t101 + t59 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, 0, t23, -t22, 0, t4 * t51 - t74 * t53 - g(3) + (-t10 * t51 + t11 * t53) * qJD(3), 0, 0, 0, 0, 0, 0, t23, 0, t22, t1 * t51 - t2 * t53 - g(3) + (t51 * t6 + t53 * t7) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, -t24, t43, t32, t83, qJDD(3), -t28 * t91 + t61, g(3) * t51 + t76 * qJD(3) + (-qJD(1) * t28 + t71) * t53 - t79, 0, 0, -t32, t43, t24, qJDD(3), -t83, t32, (2 * t89) - qJDD(4) + (t19 * t53 - t51 * t8) * qJD(1) + t61, -t67 * qJDD(1), (2 * t80) + (qJD(1) * t19 - g(3)) * t51 + (qJD(1) * t8 - t71) * t53 + (0.2e1 * qJD(4) - t76) * qJD(3) + t79, -t2 * pkin(3) - g(3) * t68 + t1 * qJ(4) + t104 * t7 - t6 * t11 - t8 * t19 + t71 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(3) - t32, t43, -t47 * t56 - t55, -t7 * qJD(3) + t8 * t91 + t103 - t64;];
tau_reg = t5;
