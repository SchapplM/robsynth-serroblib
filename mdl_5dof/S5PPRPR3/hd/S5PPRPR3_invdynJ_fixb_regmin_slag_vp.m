% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PPRPR3
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% 
% Output:
% tau_reg [5x13]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PPRPR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:05:24
% EndTime: 2019-12-05 15:05:26
% DurationCPUTime: 0.45s
% Computational Cost: add. (369->103), mult. (833->165), div. (0->0), fcn. (715->12), ass. (0->78)
t51 = cos(pkin(8));
t35 = -t51 * qJDD(1) + qJDD(4);
t48 = sin(pkin(8));
t49 = sin(pkin(7));
t52 = cos(pkin(7));
t68 = (g(1) * t52 + g(2) * t49) * t48;
t101 = -g(3) * t51 - t35 + t68;
t54 = sin(qJ(3));
t56 = cos(qJ(3));
t85 = qJD(1) * t48;
t100 = t56 * qJD(2) - t54 * t85;
t47 = sin(pkin(9));
t50 = cos(pkin(9));
t28 = t47 * t54 - t50 * t56;
t99 = g(3) * t48;
t31 = t54 * qJD(2) + t56 * t85;
t41 = t56 * qJDD(2);
t79 = t48 * qJDD(1);
t13 = qJDD(3) * pkin(3) - t31 * qJD(3) - t54 * t79 + t41;
t78 = t54 * qJDD(2);
t14 = qJD(3) * t100 + t56 * t79 + t78;
t4 = t47 * t13 + t50 * t14;
t98 = t47 * t31;
t96 = t49 * t51;
t95 = t49 * t54;
t94 = t49 * t56;
t93 = t50 * t31;
t44 = qJ(3) + pkin(9);
t39 = sin(t44);
t91 = t52 * t39;
t40 = cos(t44);
t90 = t52 * t40;
t89 = t52 * t54;
t88 = t52 * t56;
t53 = sin(qJ(5));
t55 = cos(qJ(5));
t87 = t53 * t55;
t45 = t53 ^ 2;
t86 = -t55 ^ 2 + t45;
t84 = qJD(3) * t28;
t83 = qJD(3) * t48;
t23 = t28 * t48;
t82 = qJD(5) * t23;
t81 = qJDD(1) - g(3);
t80 = t23 * qJDD(5);
t77 = t55 * qJDD(3);
t76 = qJD(3) * qJD(5);
t74 = -g(1) * t49 + g(2) * t52;
t72 = t48 * t81;
t3 = t50 * t13 - t47 * t14;
t25 = qJD(3) * pkin(3) + t100;
t7 = t50 * t25 - t98;
t29 = t47 * t56 + t50 * t54;
t58 = qJD(3) ^ 2;
t70 = t56 * qJDD(3) - t58 * t54;
t69 = -qJDD(3) * t54 - t58 * t56;
t16 = t29 * t83;
t22 = t29 * t48;
t66 = t22 * qJD(3) + qJD(5) * t51 + t16;
t26 = t29 * qJD(3);
t57 = qJD(5) ^ 2;
t65 = qJD(3) * t26 + qJDD(3) * t28 + t29 * t57;
t15 = t28 * t83;
t64 = t15 * qJD(3) - t22 * qJDD(3) - t51 * qJDD(5);
t63 = 0.2e1 * t84 * qJD(5) - qJDD(5) * t29;
t11 = t100 * t50 - t98;
t37 = t47 * pkin(3) + pkin(6);
t38 = -t50 * pkin(3) - pkin(4);
t5 = -qJD(3) * pkin(4) - t7;
t62 = -qJDD(5) * t37 + (qJD(3) * t38 + t11 + t5) * qJD(5);
t61 = -qJDD(3) * pkin(6) + g(1) * (t49 * t39 + t51 * t90) + g(2) * (t40 * t96 - t91) - t5 * qJD(3) + t40 * t99 - t4;
t60 = -g(1) * (-t51 * t89 + t94) - g(2) * (-t51 * t95 - t88);
t10 = t100 * t47 + t93;
t59 = -g(1) * (t49 * t40 - t51 * t91) - g(2) * (-t39 * t96 - t90) + qJD(3) * t10 - t37 * t57 + t39 * t99 + t3 + (pkin(4) - t38) * qJDD(3);
t34 = qJDD(5) * t55 - t57 * t53;
t33 = qJDD(5) * t53 + t57 * t55;
t8 = t47 * t25 + t93;
t1 = [t81, -g(3) + (t48 ^ 2 + t51 ^ 2) * qJDD(1), 0, t69 * t48, -t70 * t48, t7 * t15 - t8 * t16 - t3 * t22 - t4 * t23 - t35 * t51 - g(3), 0, 0, 0, 0, 0, t53 * t80 + t64 * t55 + (t53 * t66 + t55 * t82) * qJD(5), t55 * t80 - t64 * t53 + (-t53 * t82 + t55 * t66) * qJD(5); 0, qJDD(2) + t74, 0, t70, t69, -t7 * t26 - t3 * t28 + t4 * t29 - t8 * t84 + t74, 0, 0, 0, 0, 0, t53 * t63 - t55 * t65, t53 * t65 + t55 * t63; 0, 0, qJDD(3), -t54 * t72 + t41 + t60, -t78 - g(1) * (-t51 * t88 - t95) - g(2) * (-t51 * t94 + t89) - t56 * t72, t7 * t10 - t8 * t11 + (t3 * t50 + t4 * t47 + t54 * t99 + t60) * pkin(3), t45 * qJDD(3) + 0.2e1 * t76 * t87, 0.2e1 * t53 * t77 - 0.2e1 * t86 * t76, t33, t34, 0, t53 * t62 + t55 * t59, -t53 * t59 + t55 * t62; 0, 0, 0, 0, 0, -t51 * t81 + qJDD(4) - t68, 0, 0, 0, 0, 0, t34, -t33; 0, 0, 0, 0, 0, 0, -t58 * t87, t86 * t58, t53 * qJDD(3), t77, qJDD(5), -t101 * t55 + t61 * t53, t101 * t53 + t61 * t55;];
tau_reg = t1;
