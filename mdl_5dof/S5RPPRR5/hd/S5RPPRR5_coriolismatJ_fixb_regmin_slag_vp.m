% Calculate minimal parameter regressor of coriolis matrix for
% S5RPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x17]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPPRR5_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR5_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR5_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR5_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:56:40
% EndTime: 2019-12-31 17:56:41
% DurationCPUTime: 0.35s
% Computational Cost: add. (320->66), mult. (480->86), div. (0->0), fcn. (385->6), ass. (0->48)
t26 = sin(pkin(8)) * pkin(1) + qJ(3);
t32 = sin(qJ(4));
t34 = cos(qJ(4));
t40 = -cos(pkin(8)) * pkin(1) - pkin(2) - pkin(3);
t14 = t34 * t26 + t32 * t40;
t45 = qJD(1) - qJD(4);
t61 = t45 * t14;
t31 = sin(qJ(5));
t33 = cos(qJ(5));
t23 = -t31 ^ 2 + t33 ^ 2;
t60 = t45 * t23;
t13 = t32 * t26 - t34 * t40;
t11 = pkin(4) + t13;
t59 = t11 + t13;
t56 = pkin(4) * t31;
t55 = qJD(1) * t11;
t54 = qJD(1) * t33;
t53 = qJD(4) * t33;
t29 = qJD(5) * t31;
t30 = qJD(5) * t33;
t52 = qJD(5) * t34;
t51 = t23 * qJD(5);
t50 = t26 * qJD(1);
t49 = t32 * qJD(1);
t48 = t32 * qJD(3);
t47 = t34 * qJD(1);
t46 = t34 * qJD(3);
t44 = t31 * t49;
t43 = t33 * t49;
t19 = t45 * t32;
t42 = t45 * t33;
t20 = t45 * t34;
t41 = t13 / 0.2e1 - pkin(4) / 0.2e1 - t11 / 0.2e1;
t39 = -t46 + t55;
t38 = t14 * qJD(1) + t48;
t37 = t14 * qJD(4) + t48;
t1 = t41 * t31;
t36 = t1 * qJD(1) + qJD(4) * t56;
t2 = t41 * t33;
t35 = pkin(4) * t53 + t2 * qJD(1);
t25 = t31 * t30;
t17 = (-t53 + t54) * t31;
t12 = -pkin(7) + t14;
t10 = -t31 * t52 + t32 * t42;
t9 = t31 * t19 + t33 * t52;
t4 = (pkin(4) + t59) * t33 / 0.2e1;
t3 = t56 / 0.2e1 + t59 * t31 / 0.2e1;
t5 = [0, 0, 0, 0, 0, qJD(3), t26 * qJD(3), 0, t37, -t13 * qJD(4) + t46, t25, t51, 0, 0, 0, -t11 * t29 + t37 * t33, -t11 * t30 - t37 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, qJD(1), t50, 0, t49, t47, 0, 0, 0, 0, 0, t43, -t44; 0, 0, 0, 0, 0, 0, 0, 0, t61, -t45 * t13, -t25, -t51, 0, 0, 0, t3 * qJD(5) + t14 * t42, t4 * qJD(5) - t31 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, t60, -t30, t29, 0, t3 * qJD(4) - t12 * t30 - t31 * t55, t4 * qJD(4) - t11 * t54 + t12 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, t30; 0, 0, 0, 0, 0, -qJD(1), -t50, 0, -t19, -t20, 0, 0, 0, 0, 0, -t10, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t19, t20, 0, 0, 0, 0, 0, t10, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31 * t20 - t32 * t30, t33 * t20 + t32 * t29; 0, 0, 0, 0, 0, 0, 0, 0, -t38, t13 * qJD(1) - t46, -t25, -t51, 0, 0, 0, -t1 * qJD(5) - t38 * t33, -t2 * qJD(5) + t38 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t47, 0, 0, 0, 0, 0, -t43, t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t51, 0, 0, 0, -pkin(4) * t29, -pkin(4) * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, -t60, t30, -t29, 0, -pkin(7) * t30 - t36, pkin(7) * t29 - t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, -t60, 0, 0, 0, t1 * qJD(4) + t39 * t31, t2 * qJD(4) + t39 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31 * t47, -t33 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, t60, 0, 0, 0, t36, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t5;
