% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRPPR3
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRPPR3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR3_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:26:55
% EndTime: 2019-12-05 15:26:57
% DurationCPUTime: 0.33s
% Computational Cost: add. (734->92), mult. (1211->132), div. (0->0), fcn. (766->8), ass. (0->65)
t71 = pkin(3) + pkin(6);
t50 = sin(qJ(5));
t42 = t50 ^ 2;
t55 = qJD(2) ^ 2;
t70 = t42 * t55;
t52 = cos(qJ(5));
t43 = t52 ^ 2;
t69 = t43 * t55;
t62 = t52 * t55 * t50;
t33 = qJDD(5) + t62;
t68 = t50 * t33;
t34 = qJDD(5) - t62;
t67 = t52 * t34;
t47 = sin(pkin(7));
t49 = cos(pkin(7));
t32 = -t49 * g(1) - t47 * g(2);
t44 = -g(3) + qJDD(1);
t51 = sin(qJ(2));
t53 = cos(qJ(2));
t23 = -t51 * t32 + t53 * t44;
t21 = qJDD(2) * pkin(2) + t23;
t24 = t53 * t32 + t51 * t44;
t22 = -t55 * pkin(2) + t24;
t46 = sin(pkin(8));
t48 = cos(pkin(8));
t12 = t46 * t21 + t48 * t22;
t66 = t42 + t43;
t65 = t50 * qJDD(2);
t64 = qJD(2) * qJD(5);
t63 = qJDD(2) * qJ(4);
t11 = t48 * t21 - t46 * t22;
t29 = t46 * qJDD(2) + t48 * t55;
t61 = -pkin(2) * t29 - t12;
t27 = -t47 * g(1) + t49 * g(2) + qJDD(3);
t45 = qJDD(2) * pkin(3);
t10 = -t55 * qJ(4) + qJDD(4) - t11 - t45;
t56 = -qJDD(2) * pkin(6) + t10;
t5 = t50 * t27 - t52 * t56;
t6 = t52 * t27 + t50 * t56;
t2 = -t52 * t5 + t50 * t6;
t39 = t52 * qJDD(2);
t26 = -0.2e1 * t50 * t64 + t39;
t28 = t48 * qJDD(2) - t46 * t55;
t60 = t53 * t28 - t51 * t29;
t59 = t51 * t28 + t53 * t29;
t58 = pkin(2) * t28 + t11;
t40 = 2 * qJD(4) * qJD(2);
t57 = t40 + t63 + t12;
t25 = 0.2e1 * t52 * t64 + t65;
t54 = qJD(5) ^ 2;
t36 = -t54 - t69;
t35 = -t54 - t70;
t31 = t66 * t55;
t30 = t66 * qJDD(2);
t20 = t52 * t36 - t68;
t19 = t50 * t35 + t67;
t15 = t48 * t30 - t46 * t31;
t14 = -t48 * t20 + t46 * t26;
t13 = -t48 * t19 + t46 * t25;
t9 = -t55 * pkin(3) + t57;
t8 = -t71 * t55 + t57;
t4 = t48 * t11 + t46 * t12;
t3 = -t48 * t10 + t46 * t9;
t1 = -t48 * t2 + t46 * t8;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, 0, 0, 0, 0, 0, t53 * qJDD(2) - t51 * t55, -t51 * qJDD(2) - t53 * t55, 0, t53 * t23 + t51 * t24, 0, 0, 0, 0, 0, 0, t60, -t59, 0, t51 * (-t46 * t11 + t48 * t12) + t53 * t4, 0, 0, 0, 0, 0, 0, 0, -t60, t59, t51 * (t46 * t10 + t48 * t9) + t53 * t3, 0, 0, 0, 0, 0, 0, t51 * (t46 * t19 + t48 * t25) + t53 * t13, t51 * (t46 * t20 + t48 * t26) + t53 * t14, t51 * (-t46 * t30 - t48 * t31) + t53 * t15, t51 * (t46 * t2 + t48 * t8) + t53 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t23, -t24, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t58, t61, 0, pkin(2) * t4, qJDD(2), 0, 0, 0, 0, 0, 0, qJDD(4) - 0.2e1 * t45 - t58, t40 - t61 + 0.2e1 * t63, pkin(2) * t3 - pkin(3) * t10 + qJ(4) * t9, t26 * t52, -t52 * t25 - t50 * t26, t67 - t50 * (t54 - t69), t25 * t50, t52 * (-t54 + t70) - t68, 0, pkin(2) * t13 + qJ(4) * t25 - t71 * t19 + t50 * t8, pkin(2) * t14 + qJ(4) * t26 - t71 * t20 + t52 * t8, pkin(2) * t15 - qJ(4) * t31 + t71 * t30 - t2, pkin(2) * t1 + qJ(4) * t8 - t71 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, 0, 0, 0, 0, 0, -t50 * t34 + t52 * t35, -t52 * t33 - t50 * t36, 0, t50 * t5 + t52 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t55, t10, 0, 0, 0, 0, 0, 0, t19, t20, -t30, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, (-t42 + t43) * t55, t39, -t62, -t65, qJDD(5), -t5, -t6, 0, 0;];
tauJ_reg = t7;
