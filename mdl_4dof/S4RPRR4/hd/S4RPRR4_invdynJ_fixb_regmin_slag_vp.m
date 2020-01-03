% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RPRR4
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
% tau_reg [4x18]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPRR4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR4_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:50:37
% EndTime: 2019-12-31 16:50:39
% DurationCPUTime: 0.64s
% Computational Cost: add. (560->149), mult. (1199->233), div. (0->0), fcn. (781->10), ass. (0->89)
t50 = sin(qJ(3));
t91 = qJD(4) * t50;
t108 = -qJD(1) * t91 + qJDD(3);
t49 = sin(qJ(4));
t52 = cos(qJ(4));
t83 = t50 * qJDD(1);
t53 = cos(qJ(3));
t85 = t53 * qJD(1);
t7 = ((qJD(4) + t85) * qJD(3) + t83) * t49 - t108 * t52;
t34 = -qJD(4) + t85;
t86 = t52 * qJD(3);
t79 = t53 * t86;
t63 = -t49 * t91 + t79;
t42 = t53 * qJDD(1);
t82 = qJD(1) * qJD(3);
t75 = t50 * t82;
t20 = qJDD(4) - t42 + t75;
t98 = t52 * t20;
t107 = -t63 * t34 + t50 * t98;
t47 = sin(pkin(7));
t35 = t47 * pkin(1) + pkin(5);
t28 = t35 * qJDD(1);
t30 = t35 * qJD(1);
t17 = t50 * qJD(2) + t53 * t30;
t88 = t17 * qJD(3);
t4 = -qJDD(3) * pkin(3) - t53 * qJDD(2) + t50 * t28 + t88;
t44 = qJ(1) + pkin(7);
t38 = sin(t44);
t39 = cos(t44);
t69 = g(1) * t39 + g(2) * t38;
t70 = pkin(3) * t50 - pkin(6) * t53;
t105 = (pkin(6) * qJD(4) + t70 * qJD(1)) * t34 + t69 * t50 - g(3) * t53 - t4;
t6 = qJD(1) * t79 + qJD(4) * t86 + t108 * t49 + t52 * t83;
t104 = t6 * t49;
t95 = qJD(1) * t50;
t21 = t49 * t95 - t86;
t103 = t21 * t34;
t87 = t49 * qJD(3);
t23 = t52 * t95 + t87;
t102 = t23 * t34;
t101 = t34 * t52;
t100 = t49 * t20;
t99 = t49 * t53;
t97 = t52 * t53;
t45 = t50 ^ 2;
t96 = -t53 ^ 2 + t45;
t48 = cos(pkin(7));
t36 = -t48 * pkin(1) - pkin(2);
t31 = qJD(1) * t36;
t94 = qJD(3) * t21;
t93 = qJD(3) * t50;
t92 = qJD(4) * t49;
t90 = qJD(4) * t52;
t16 = t53 * qJD(2) - t50 * t30;
t89 = t16 * qJD(3);
t84 = qJDD(2) - g(3);
t80 = t34 * t87;
t77 = t23 * t93 - t6 * t53;
t10 = qJD(3) * pkin(6) + t17;
t74 = t34 * t35 + t10;
t19 = -t53 * pkin(3) - t50 * pkin(6) + t36;
t11 = t19 * qJD(1);
t72 = -qJDD(3) * pkin(6) - qJD(4) * t11 - t50 * qJDD(2) - t53 * t28 - t89;
t51 = sin(qJ(1));
t54 = cos(qJ(1));
t68 = g(1) * t51 - g(2) * t54;
t67 = qJD(3) * t30 - t84;
t65 = g(3) * t50 + t72;
t64 = t34 * t90 - t100;
t25 = t70 * qJD(3);
t9 = -qJD(3) * pkin(3) - t16;
t62 = qJD(3) * t9 - t20 * t35 - t72;
t61 = 0.2e1 * t31 * qJD(3) - qJDD(3) * t35;
t60 = -pkin(6) * t20 + (-t16 - t9) * t34;
t55 = qJD(3) ^ 2;
t59 = g(1) * t38 - g(2) * t39 - 0.2e1 * qJDD(1) * t36 - t35 * t55;
t58 = -qJD(1) * t31 - qJD(2) * qJD(3) - t28 + t69;
t56 = qJD(1) ^ 2;
t27 = qJDD(3) * t53 - t55 * t50;
t26 = qJDD(3) * t50 + t55 * t53;
t15 = t38 * t49 + t39 * t97;
t14 = t38 * t52 - t39 * t99;
t13 = -t38 * t97 + t39 * t49;
t12 = t38 * t99 + t39 * t52;
t8 = qJD(1) * t25 + t19 * qJDD(1);
t5 = t52 * t8;
t2 = t52 * t10 + t49 * t11;
t1 = -t49 * t10 + t52 * t11;
t3 = [qJDD(1), t68, g(1) * t54 + g(2) * t51, (t68 + (t47 ^ 2 + t48 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t45 * qJDD(1) + 0.2e1 * t53 * t75, 0.2e1 * t50 * t42 - 0.2e1 * t96 * t82, t26, t27, 0, t61 * t50 + t59 * t53, -t59 * t50 + t61 * t53, t6 * t52 * t50 + t63 * t23, (-t21 * t52 - t23 * t49) * t53 * qJD(3) + (-t104 - t52 * t7 + (t21 * t49 - t23 * t52) * qJD(4)) * t50, t107 + t77, (t7 + t80) * t53 + (t64 - t94) * t50, -t20 * t53 - t34 * t93, t50 * t35 * t7 - g(1) * t13 - g(2) * t15 - t5 * t53 + (t53 * t35 * t21 + t1 * t50) * qJD(3) + (t19 * t20 - t25 * t34 + (t9 * t50 + t74 * t53) * qJD(4)) * t52 + (-(-qJD(4) * t19 + t35 * t93) * t34 + t4 * t50 + t62 * t53) * t49, (t19 * t90 + t49 * t25) * t34 - t19 * t100 - g(1) * t12 - g(2) * t14 + (qJD(3) * t35 * t23 + (-t74 * qJD(4) + t8) * t49 + t62 * t52) * t53 + (-t9 * t92 + t35 * t6 + t4 * t52 + (-t35 * t101 - t2) * qJD(3)) * t50; 0, 0, 0, t84, 0, 0, 0, 0, 0, t27, -t26, 0, 0, 0, 0, 0, (-t7 + t80) * t53 + (t64 + t94) * t50, -t107 + t77; 0, 0, 0, 0, -t50 * t56 * t53, t96 * t56, t83, t42, qJDD(3), t58 * t50 - t67 * t53 + t88, t67 * t50 + t58 * t53 + t89, -t23 * t101 + t104, (t6 + t103) * t52 + (-t7 + t102) * t49, (-t23 * t50 + t34 * t97) * qJD(1) - t64, t34 * t92 + t98 + (t21 * t50 - t34 * t99) * qJD(1), t34 * t95, -pkin(3) * t7 - t1 * t95 + t105 * t52 - t17 * t21 + t60 * t49, -pkin(3) * t6 - t105 * t49 - t17 * t23 + t2 * t95 + t60 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23 * t21, -t21 ^ 2 + t23 ^ 2, t6 - t103, -t102 - t7, t20, -g(1) * t14 + g(2) * t12 - t10 * t90 - t2 * t34 - t9 * t23 + t65 * t49 + t5, g(1) * t15 - g(2) * t13 - t1 * t34 + t9 * t21 + (qJD(4) * t10 - t8) * t49 + t65 * t52;];
tau_reg = t3;
