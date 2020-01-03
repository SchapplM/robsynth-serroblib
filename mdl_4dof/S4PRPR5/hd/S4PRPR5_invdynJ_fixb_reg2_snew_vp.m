% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4PRPR5
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
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4PRPR5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR5_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR5_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR5_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:23:18
% EndTime: 2019-12-31 16:23:19
% DurationCPUTime: 0.21s
% Computational Cost: add. (460->75), mult. (815->116), div. (0->0), fcn. (546->8), ass. (0->52)
t48 = sin(qJ(4));
t50 = cos(qJ(4));
t53 = qJD(2) ^ 2;
t34 = t48 * t53 * t50;
t30 = qJDD(4) + t34;
t59 = t48 * t30;
t31 = qJDD(4) - t34;
t58 = t50 * t31;
t44 = sin(pkin(6));
t46 = cos(pkin(6));
t29 = -t46 * g(1) - t44 * g(2);
t41 = -g(3) + qJDD(1);
t49 = sin(qJ(2));
t51 = cos(qJ(2));
t21 = t51 * t29 + t49 * t41;
t19 = -t53 * pkin(2) + t21;
t43 = sin(pkin(7));
t45 = cos(pkin(7));
t20 = -t49 * t29 + t51 * t41;
t54 = qJDD(2) * pkin(2) + t20;
t10 = t45 * t19 + t43 * t54;
t57 = t48 * qJDD(2);
t56 = qJD(2) * qJD(4);
t55 = -t44 * g(1) + t46 * g(2) + qJDD(3);
t8 = -t53 * pkin(3) + qJDD(2) * pkin(5) + t10;
t5 = t48 * t8 - t50 * t55;
t6 = t48 * t55 + t50 * t8;
t2 = t48 * t5 + t50 * t6;
t9 = -t43 * t19 + t45 * t54;
t36 = t50 * qJDD(2);
t24 = -0.2e1 * t48 * t56 + t36;
t23 = 0.2e1 * t50 * t56 + t57;
t52 = qJD(4) ^ 2;
t40 = t50 ^ 2;
t39 = t48 ^ 2;
t38 = t40 * t53;
t37 = t39 * t53;
t33 = -t38 - t52;
t32 = -t37 - t52;
t28 = t37 + t38;
t27 = (t39 + t40) * qJDD(2);
t26 = -t43 * qJDD(2) - t45 * t53;
t25 = t45 * qJDD(2) - t43 * t53;
t18 = -t48 * t32 - t58;
t17 = t50 * t33 - t59;
t13 = t43 * t27 + t45 * t28;
t12 = t43 * t18 - t45 * t23;
t11 = t43 * t17 + t45 * t24;
t7 = -qJDD(2) * pkin(3) - t53 * pkin(5) - t9;
t3 = t43 * t10 + t45 * t9;
t1 = t43 * t2 - t45 * t7;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, 0, 0, 0, 0, 0, t51 * qJDD(2) - t49 * t53, -t49 * qJDD(2) - t51 * t53, 0, t51 * t20 + t49 * t21, 0, 0, 0, 0, 0, 0, t51 * t25 + t49 * t26, -t49 * t25 + t51 * t26, 0, t49 * (t45 * t10 - t43 * t9) + t51 * t3, 0, 0, 0, 0, 0, 0, t49 * (t45 * t17 - t43 * t24) + t51 * t11, t49 * (t45 * t18 + t43 * t23) + t51 * t12, t49 * (t45 * t27 - t43 * t28) + t51 * t13, t49 * (t45 * t2 + t43 * t7) + t51 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t20, -t21, 0, 0, 0, 0, 0, 0, 0, qJDD(2), pkin(2) * t25 + t9, pkin(2) * t26 - t10, 0, pkin(2) * t3, t23 * t48, t50 * t23 + t48 * t24, t59 + t50 * (-t37 + t52), t24 * t50, t48 * (t38 - t52) + t58, 0, pkin(2) * t11 + pkin(3) * t24 + pkin(5) * t17 - t50 * t7, pkin(2) * t12 - pkin(3) * t23 + pkin(5) * t18 + t48 * t7, pkin(2) * t13 + pkin(3) * t28 + pkin(5) * t27 + t2, pkin(2) * t1 - pkin(3) * t7 + pkin(5) * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, 0, 0, 0, 0, 0, t50 * t30 + t48 * t33, -t48 * t31 + t50 * t32, 0, t48 * t6 - t50 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, t37 - t38, t57, t34, t36, qJDD(4), -t5, -t6, 0, 0;];
tauJ_reg = t4;
