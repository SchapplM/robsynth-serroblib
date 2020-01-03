% Calculate minimal parameter regressor of coriolis matrix for
% S5PRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x17]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRPRR9_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR9_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR9_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR9_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:39:46
% EndTime: 2019-12-31 17:39:47
% DurationCPUTime: 0.33s
% Computational Cost: add. (222->64), mult. (381->84), div. (0->0), fcn. (286->4), ass. (0->47)
t31 = sin(qJ(4));
t33 = cos(qJ(4));
t55 = -pkin(2) - pkin(3);
t17 = t33 * qJ(3) + t31 * t55;
t43 = qJD(2) - qJD(4);
t60 = t43 * t17;
t30 = sin(qJ(5));
t32 = cos(qJ(5));
t22 = -t30 ^ 2 + t32 ^ 2;
t59 = t43 * t22;
t16 = t31 * qJ(3) - t33 * t55;
t14 = pkin(4) + t16;
t58 = t14 + t16;
t54 = pkin(4) * t30;
t53 = qJD(2) * t14;
t52 = qJD(2) * t32;
t51 = qJD(4) * t32;
t28 = qJD(5) * t30;
t29 = qJD(5) * t32;
t50 = qJD(5) * t33;
t49 = t22 * qJD(5);
t48 = t31 * qJD(2);
t47 = t31 * qJD(3);
t46 = t33 * qJD(2);
t45 = t33 * qJD(3);
t44 = qJ(3) * qJD(2);
t42 = t30 * t48;
t41 = t32 * t48;
t18 = t43 * t31;
t40 = t43 * t32;
t19 = t43 * t33;
t39 = t16 / 0.2e1 - pkin(4) / 0.2e1 - t14 / 0.2e1;
t38 = -t45 + t53;
t37 = t17 * qJD(2) + t47;
t36 = t17 * qJD(4) + t47;
t1 = t39 * t30;
t35 = t1 * qJD(2) + qJD(4) * t54;
t2 = t39 * t32;
t34 = pkin(4) * t51 + t2 * qJD(2);
t24 = t30 * t29;
t15 = -pkin(7) + t17;
t13 = (-t51 + t52) * t30;
t6 = -t30 * t50 + t31 * t40;
t5 = t30 * t18 + t32 * t50;
t4 = (pkin(4) + t58) * t32 / 0.2e1;
t3 = t54 / 0.2e1 + t58 * t30 / 0.2e1;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, qJD(3), qJ(3) * qJD(3), 0, t36, -t16 * qJD(4) + t45, t24, t49, 0, 0, 0, -t14 * t28 + t36 * t32, -t14 * t29 - t36 * t30; 0, 0, 0, 0, 0, qJD(2), t44, 0, t48, t46, 0, 0, 0, 0, 0, t41, -t42; 0, 0, 0, 0, 0, 0, 0, 0, t60, -t43 * t16, -t24, -t49, 0, 0, 0, t3 * qJD(5) + t17 * t40, t4 * qJD(5) - t30 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t59, -t29, t28, 0, t3 * qJD(4) - t15 * t29 - t30 * t53, t4 * qJD(4) - t14 * t52 + t15 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -qJD(2), -t44, 0, -t18, -t19, 0, 0, 0, 0, 0, -t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t18, t19, 0, 0, 0, 0, 0, t6, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30 * t19 - t31 * t29, t32 * t19 + t31 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t37, t16 * qJD(2) - t45, -t24, -t49, 0, 0, 0, -t1 * qJD(5) - t37 * t32, -t2 * qJD(5) + t37 * t30; 0, 0, 0, 0, 0, 0, 0, 0, -t48, -t46, 0, 0, 0, 0, 0, -t41, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t49, 0, 0, 0, -pkin(4) * t28, -pkin(4) * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t59, t29, -t28, 0, -pkin(7) * t29 - t35, pkin(7) * t28 - t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t59, 0, 0, 0, t1 * qJD(4) + t38 * t30, t2 * qJD(4) + t38 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30 * t46, -t32 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t59, 0, 0, 0, t35, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t7;
