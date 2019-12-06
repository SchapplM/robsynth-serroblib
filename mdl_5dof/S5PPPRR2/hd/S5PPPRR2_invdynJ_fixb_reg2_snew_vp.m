% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PPPRR2
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
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PPPRR2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR2_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPPRR2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:59:37
% EndTime: 2019-12-05 14:59:39
% DurationCPUTime: 0.31s
% Computational Cost: add. (753->82), mult. (1217->131), div. (0->0), fcn. (998->10), ass. (0->66)
t55 = sin(pkin(7));
t58 = cos(pkin(7));
t38 = -t58 * g(1) - t55 * g(2);
t51 = -g(3) + qJDD(1);
t54 = sin(pkin(8));
t57 = cos(pkin(8));
t30 = t57 * t38 + t54 * t51;
t53 = sin(pkin(9));
t56 = cos(pkin(9));
t66 = -t55 * g(1) + t58 * g(2) + qJDD(2);
t20 = t53 * t30 - t56 * t66;
t73 = t53 * t20;
t72 = t54 * t56;
t60 = sin(qJ(5));
t62 = cos(qJ(5));
t65 = qJD(4) ^ 2;
t43 = t60 * t65 * t62;
t39 = qJDD(5) + t43;
t71 = t60 * t39;
t40 = qJDD(5) - t43;
t70 = t62 * t40;
t69 = t60 * qJDD(4);
t68 = qJD(4) * qJD(5);
t21 = t56 * t30 + t53 * t66;
t67 = -t54 * t38 + t57 * t51;
t29 = qJDD(3) - t67;
t61 = sin(qJ(4));
t63 = cos(qJ(4));
t13 = t63 * t21 + t61 * t29;
t11 = -t65 * pkin(4) + qJDD(4) * pkin(6) + t13;
t8 = t60 * t11 - t62 * t20;
t9 = t62 * t11 + t60 * t20;
t4 = t60 * t8 + t62 * t9;
t12 = -t61 * t21 + t63 * t29;
t46 = t62 * qJDD(4);
t33 = -0.2e1 * t60 * t68 + t46;
t32 = 0.2e1 * t62 * t68 + t69;
t64 = qJD(5) ^ 2;
t50 = t62 ^ 2;
t49 = t60 ^ 2;
t48 = t50 * t65;
t47 = t49 * t65;
t42 = -t48 - t64;
t41 = -t47 - t64;
t37 = t47 + t48;
t36 = -t63 * qJDD(4) + t61 * t65;
t35 = t61 * qJDD(4) + t63 * t65;
t34 = (t49 + t50) * qJDD(4);
t28 = -t60 * t41 - t70;
t27 = t62 * t42 - t71;
t26 = t60 * t40 - t62 * t41;
t25 = -t62 * t39 - t60 * t42;
t23 = t63 * t34 - t61 * t37;
t22 = t61 * t34 + t63 * t37;
t19 = t63 * t28 + t61 * t32;
t18 = t63 * t27 - t61 * t33;
t17 = t61 * t28 - t63 * t32;
t16 = t61 * t27 + t63 * t33;
t14 = t56 * t20;
t10 = -qJDD(4) * pkin(4) - t65 * pkin(6) - t12;
t6 = -t61 * t12 + t63 * t13;
t5 = t63 * t12 + t61 * t13;
t3 = -t60 * t9 + t62 * t8;
t2 = t61 * t10 + t63 * t4;
t1 = -t63 * t10 + t61 * t4;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54 * t30 + t57 * t67, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54 * (t56 * t21 + t73) - t57 * t29, 0, 0, 0, 0, 0, 0, -t35 * t72 + t57 * t36, t57 * t35 + t36 * t72, 0, t54 * (t56 * t6 + t73) - t57 * t5, 0, 0, 0, 0, 0, 0, t54 * (t56 * t18 - t53 * t25) - t57 * t16, t54 * (t56 * t19 - t53 * t26) - t57 * t17, -t57 * t22 + t23 * t72, t54 * (t56 * t2 - t53 * t3) - t57 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53 * t21 - t14, 0, 0, 0, 0, 0, 0, -t53 * t35, t53 * t36, 0, t53 * t6 - t14, 0, 0, 0, 0, 0, 0, t53 * t18 + t56 * t25, t53 * t19 + t56 * t26, t53 * t23, t53 * t2 + t56 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, 0, 0, 0, 0, 0, -t36, -t35, 0, t5, 0, 0, 0, 0, 0, 0, t16, t17, t22, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(4), t12, -t13, 0, 0, t32 * t60, t62 * t32 + t60 * t33, t71 + t62 * (-t47 + t64), t33 * t62, t60 * (t48 - t64) + t70, 0, pkin(4) * t33 + pkin(6) * t27 - t62 * t10, -pkin(4) * t32 + pkin(6) * t28 + t60 * t10, pkin(4) * t37 + pkin(6) * t34 + t4, -pkin(4) * t10 + pkin(6) * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, t47 - t48, t69, t43, t46, qJDD(5), -t8, -t9, 0, 0;];
tauJ_reg = t7;
