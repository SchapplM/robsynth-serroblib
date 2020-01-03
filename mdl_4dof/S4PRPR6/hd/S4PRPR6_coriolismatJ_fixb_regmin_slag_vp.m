% Calculate minimal parameter regressor of coriolis matrix for
% S4PRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x15]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4PRPR6_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR6_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:24:42
% EndTime: 2019-12-31 16:24:43
% DurationCPUTime: 0.22s
% Computational Cost: add. (134->50), mult. (384->85), div. (0->0), fcn. (377->6), ass. (0->45)
t30 = sin(pkin(7));
t28 = t30 ^ 2;
t31 = cos(pkin(7));
t29 = t31 ^ 2;
t22 = t29 + t28;
t55 = cos(qJ(4));
t32 = sin(qJ(4));
t54 = t32 * t30;
t53 = t32 * t31;
t52 = pkin(5) + qJ(3);
t39 = t55 * t31;
t12 = -t39 + t54;
t40 = t55 * t30;
t14 = t40 + t53;
t1 = t12 ^ 2 - t14 ^ 2;
t51 = t1 * qJD(2);
t34 = cos(qJ(2));
t36 = -t53 / 0.2e1 - t40 / 0.2e1;
t2 = (t14 / 0.2e1 + t36) * t34;
t50 = t2 * qJD(1);
t35 = -t39 / 0.2e1 + t54 / 0.2e1;
t3 = (-t12 / 0.2e1 + t35) * t34;
t49 = t3 * qJD(1);
t33 = sin(qJ(2));
t6 = (-0.1e1 + t22) * t34 * t33;
t48 = t6 * qJD(1);
t47 = t12 * qJD(2);
t11 = t12 * qJD(4);
t46 = t14 * qJD(2);
t45 = t14 * qJD(4);
t44 = t22 * qJD(2);
t43 = t33 * qJD(2);
t42 = t34 * qJD(2);
t41 = t12 * t46;
t15 = t22 * qJ(3);
t26 = -t31 * pkin(3) - pkin(2);
t38 = qJD(2) * t26 + qJD(3);
t7 = (0.1e1 / 0.2e1 - t29 / 0.2e1 - t28 / 0.2e1) * t33;
t37 = t7 * qJD(1) - t15 * qJD(2);
t21 = t52 * t31;
t20 = t52 * t30;
t8 = (0.1e1 + t22) * t33 / 0.2e1;
t5 = (-t14 / 0.2e1 + t36) * t34;
t4 = (t12 / 0.2e1 + t35) * t34;
t9 = [0, 0, 0, 0, 0, 0, 0, t6 * qJD(2), 0, 0, 0, 0, 0, 0, 0; 0, 0, -t43, -t42, -t31 * t43, t30 * t43, t22 * t42, t48 + (-t33 * pkin(2) + t15 * t34) * qJD(2) + t8 * qJD(3), 0, 0, 0, 0, 0, t5 * qJD(4) + t12 * t43, t4 * qJD(4) + t14 * t43; 0, 0, 0, 0, 0, 0, 0, t8 * qJD(2), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 * qJD(2) + t11 * t33, t4 * qJD(2) + t33 * t45; 0, 0, 0, 0, 0, 0, 0, -t7 * qJD(3) - t48, 0, 0, 0, 0, 0, -t2 * qJD(4), -t3 * qJD(4); 0, 0, 0, 0, 0, 0, t22 * qJD(3), t15 * qJD(3), -t12 * t45, t1 * qJD(4), 0, 0, 0, t26 * t45, -t26 * t11; 0, 0, 0, 0, 0, 0, t44, -t37, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t41, t51, -t11, -t45, 0, -t50 + t26 * t46 + (t32 * t20 - t55 * t21) * qJD(4), -t49 - t26 * t47 + (t55 * t20 + t32 * t21) * qJD(4); 0, 0, 0, 0, 0, 0, 0, t7 * qJD(2), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -t44, t37, 0, 0, 0, 0, 0, t45, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2 * qJD(2), t3 * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, t41, -t51, 0, 0, 0, -t14 * t38 + t50, t12 * t38 + t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t9;
