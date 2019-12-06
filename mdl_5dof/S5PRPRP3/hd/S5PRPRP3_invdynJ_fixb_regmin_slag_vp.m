% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRPRP3
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% tau_reg [5x14]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRPRP3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP3_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:33:44
% EndTime: 2019-12-05 15:33:45
% DurationCPUTime: 0.57s
% Computational Cost: add. (616->133), mult. (1258->183), div. (0->0), fcn. (877->10), ass. (0->83)
t48 = qJ(2) + pkin(8);
t42 = sin(t48);
t43 = cos(t48);
t52 = sin(pkin(7));
t54 = cos(pkin(7));
t76 = g(1) * t54 + g(2) * t52;
t107 = -g(3) * t43 + t76 * t42;
t106 = qJ(5) + pkin(6);
t53 = cos(pkin(8));
t57 = sin(qJ(2));
t90 = qJD(1) * t57;
t37 = t53 * t90;
t51 = sin(pkin(8));
t59 = cos(qJ(2));
t88 = t59 * qJD(1);
t20 = t51 * t88 + t37;
t39 = t51 * pkin(2) + pkin(6);
t95 = t53 * pkin(2);
t40 = -pkin(3) - t95;
t60 = qJD(4) ^ 2;
t83 = qJD(1) * qJD(2);
t104 = t57 * qJDD(1) + t59 * t83;
t45 = t59 * qJDD(1);
t23 = qJDD(2) * pkin(2) - t57 * t83 + t45;
t9 = -t104 * t51 + t53 * t23;
t105 = qJD(2) * t20 - t39 * t60 + t9 + (pkin(3) - t40) * qJDD(2) + t107;
t102 = -g(1) * t52 + g(2) * t54;
t103 = qJDD(3) + t102;
t68 = g(3) * t42 + t76 * t43;
t56 = sin(qJ(4));
t58 = cos(qJ(4));
t33 = qJD(2) * pkin(2) + t88;
t17 = t51 * t33 + t37;
t77 = t106 * qJD(2) + t17;
t7 = t58 * qJD(3) - t77 * t56;
t92 = qJD(4) * pkin(4);
t4 = t7 + t92;
t101 = -t7 + t4;
t49 = t56 ^ 2;
t50 = t58 ^ 2;
t94 = t49 - t50;
t93 = t49 + t50;
t91 = qJ(5) + t39;
t26 = t51 * t57 - t53 * t59;
t89 = qJD(2) * t26;
t87 = qJDD(1) - g(3);
t86 = t56 * qJDD(2);
t84 = t58 * qJDD(2);
t82 = qJD(2) * qJD(4);
t10 = t104 * t53 + t51 * t23;
t41 = t58 * pkin(4) + pkin(3);
t79 = t56 * t82;
t36 = t51 * t90;
t16 = t53 * t33 - t36;
t78 = qJD(4) * t91;
t8 = t56 * qJD(3) + t77 * t58;
t74 = t4 * t56 - t8 * t58;
t27 = t51 * t59 + t53 * t57;
t73 = t77 * qJD(4);
t71 = t102 * t58;
t19 = t27 * qJD(2);
t69 = qJD(2) * t19 + qJDD(2) * t26 + t27 * t60;
t67 = 0.2e1 * t89 * qJD(4) - qJDD(4) * t27;
t66 = -g(3) * t59 + t76 * t57;
t12 = -qJD(2) * pkin(3) - t16;
t22 = t53 * t88 - t36;
t65 = -qJDD(4) * t39 + (qJD(2) * t40 + t12 + t22) * qJD(4);
t6 = qJDD(2) * pkin(6) + t10;
t64 = qJ(5) * qJDD(2) + qJD(2) * qJD(5) + qJD(4) * qJD(3) + t6;
t3 = pkin(4) * t79 - t41 * qJDD(2) + qJDD(5) - t9;
t62 = -t12 * qJD(2) - t6 + t68;
t61 = qJD(2) ^ 2;
t44 = t58 * qJDD(3);
t30 = qJDD(4) * t58 - t60 * t56;
t29 = qJDD(4) * t56 + t60 * t58;
t25 = t91 * t58;
t24 = t91 * t56;
t15 = -t56 * qJD(5) - t58 * t78;
t14 = t58 * qJD(5) - t56 * t78;
t11 = -t41 * qJD(2) + qJD(5) - t16;
t2 = (qJDD(3) - t73) * t56 + t64 * t58;
t1 = qJDD(4) * pkin(4) - t64 * t56 - t58 * t73 + t44;
t5 = [t87, 0, t59 * qJDD(2) - t61 * t57, -qJDD(2) * t57 - t61 * t59, t10 * t27 - t16 * t19 - t17 * t89 - t9 * t26 - g(3), 0, 0, 0, 0, 0, t67 * t56 - t69 * t58, t69 * t56 + t67 * t58, (-qJD(2) * t89 + qJDD(2) * t27) * t93, t11 * t19 + t3 * t26 - g(3) + t74 * t89 + (-t1 * t56 + t2 * t58 + (-t4 * t58 - t56 * t8) * qJD(4)) * t27; 0, qJDD(2), t45 + t66, -t87 * t57 + t76 * t59, t16 * t20 - t17 * t22 + (t10 * t51 + t53 * t9 + t66) * pkin(2), t49 * qJDD(2) + 0.2e1 * t58 * t79, 0.2e1 * t56 * t84 - 0.2e1 * t94 * t82, t29, t30, 0, t105 * t58 + t65 * t56, -t105 * t56 + t65 * t58, (-qJD(4) * t4 + qJDD(2) * t25 + t2) * t58 + (-qJD(4) * t8 + qJDD(2) * t24 - t1) * t56 + (t14 * t58 - t15 * t56 - t93 * t22 + (t24 * t58 - t25 * t56) * qJD(4)) * qJD(2) - t68, t2 * t25 + t8 * t14 - t1 * t24 + t4 * t15 + t3 * (-t41 - t95) - g(3) * (t59 * pkin(2) + t106 * t42 + t43 * t41) + t74 * t22 + (t56 * t92 - t20) * t11 + t76 * (pkin(2) * t57 - t106 * t43 + t41 * t42); 0, 0, 0, 0, t103, 0, 0, 0, 0, 0, t30, -t29, 0, -t74 * qJD(4) + t1 * t58 + t2 * t56 + t102; 0, 0, 0, 0, 0, -t56 * t61 * t58, t94 * t61, t86, t84, qJDD(4), t62 * t56 + t44 + t71, -t103 * t56 + t62 * t58, -pkin(4) * t86 + (-t92 + t101) * t58 * qJD(2), t101 * t8 + (t1 + t71 + (-t11 * qJD(2) + t68) * t56) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93 * t61, t74 * qJD(2) - t107 + t3;];
tau_reg = t5;
