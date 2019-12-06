% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% 
% Output:
% tauc_reg [5x13]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PPRPR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:05:22
% EndTime: 2019-12-05 15:05:23
% DurationCPUTime: 0.21s
% Computational Cost: add. (202->51), mult. (565->95), div. (0->0), fcn. (444->8), ass. (0->48)
t36 = sin(qJ(3));
t38 = cos(qJ(3));
t32 = sin(pkin(8));
t52 = qJD(1) * t32;
t64 = t38 * qJD(2) - t36 * t52;
t31 = sin(pkin(9));
t33 = cos(pkin(9));
t20 = t31 * t36 - t33 * t38;
t23 = t36 * qJD(2) + t38 * t52;
t63 = t31 * t23;
t61 = t33 * t23;
t35 = sin(qJ(5));
t37 = cos(qJ(5));
t59 = t35 * t37;
t39 = qJD(5) ^ 2;
t58 = t39 * t35;
t57 = t39 * t37;
t40 = qJD(3) ^ 2;
t56 = t40 * t36;
t55 = t40 * t38;
t54 = t35 ^ 2 - t37 ^ 2;
t51 = qJD(3) * t32;
t9 = t20 * t51;
t53 = t9 * qJD(3);
t50 = t35 * qJD(5);
t49 = t37 * qJD(5);
t47 = 0.2e1 * qJD(3) * qJD(5);
t16 = qJD(3) * pkin(3) + t64;
t5 = t33 * t16 - t63;
t1 = -qJD(3) * pkin(4) - t5;
t19 = t64 * qJD(3);
t41 = t23 * qJD(3);
t4 = t33 * t19 - t31 * t41;
t46 = -t1 * qJD(3) - t4;
t21 = t31 * t38 + t33 * t36;
t3 = t31 * t19 + t33 * t41;
t7 = t31 * t64 + t61;
t45 = qJD(3) * t7 - (t31 * pkin(3) + pkin(6)) * t39 - t3;
t8 = t33 * t64 - t63;
t44 = qJD(5) * (qJD(3) * (-t33 * pkin(3) - pkin(4)) + t1 + t8);
t10 = t21 * t51;
t13 = t21 * t32;
t42 = t13 * qJD(3) + qJD(5) * cos(pkin(8)) + t10;
t18 = t20 * qJD(3);
t17 = t21 * qJD(3);
t14 = t20 * t32;
t6 = t31 * t16 + t61;
t2 = [0, 0, 0, -t32 * t55, t32 * t56, -t6 * t10 + t3 * t13 - t4 * t14 + t5 * t9, 0, 0, 0, 0, 0, t37 * t53 + (t14 * t49 + t35 * t42) * qJD(5), -t35 * t53 + (-t14 * t50 + t37 * t42) * qJD(5); 0, 0, 0, -t56, -t55, -t5 * t17 - t6 * t18 + t3 * t20 + t4 * t21, 0, 0, 0, 0, 0, t18 * t50 - t21 * t57 + (-t17 * t37 + t20 * t50) * qJD(3), t18 * t49 + t21 * t58 + (t17 * t35 + t20 * t49) * qJD(3); 0, 0, 0, 0, 0, t5 * t7 - t6 * t8 + (-t3 * t33 + t31 * t4) * pkin(3), t47 * t59, -t54 * t47, t57, -t58, 0, t35 * t44 + t37 * t45, -t35 * t45 + t37 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, -t57; 0, 0, 0, 0, 0, 0, -t40 * t59, t54 * t40, 0, 0, 0, t46 * t35, t46 * t37;];
tauc_reg = t2;
