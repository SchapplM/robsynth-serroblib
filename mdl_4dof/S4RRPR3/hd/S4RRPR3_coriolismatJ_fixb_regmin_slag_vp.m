% Calculate minimal parameter regressor of coriolis matrix for
% S4RRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x14]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRPR3_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR3_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:01:35
% EndTime: 2019-12-31 17:01:35
% DurationCPUTime: 0.23s
% Computational Cost: add. (197->47), mult. (464->77), div. (0->0), fcn. (360->6), ass. (0->50)
t33 = cos(qJ(2));
t26 = t33 * pkin(1) + pkin(2);
t29 = cos(pkin(7));
t28 = sin(pkin(7));
t31 = sin(qJ(2));
t53 = t28 * t31;
t36 = -pkin(1) * t53 + t29 * t26;
t13 = -pkin(3) - t36;
t25 = -t29 * pkin(2) - pkin(3);
t61 = t13 / 0.2e1 + t25 / 0.2e1;
t17 = (t29 * t33 - t53) * pkin(1);
t56 = -t17 / 0.2e1;
t60 = t56 + t61;
t32 = cos(qJ(4));
t43 = qJD(1) + qJD(2);
t59 = t32 * t43;
t30 = sin(qJ(4));
t22 = -t30 ^ 2 + t32 ^ 2;
t58 = t43 * t22;
t52 = t29 * t31;
t51 = pkin(1) * qJD(1);
t50 = pkin(1) * qJD(2);
t16 = (t28 * t33 + t52) * pkin(1);
t38 = pkin(1) * t52 + t28 * t26;
t1 = -t36 * t16 + t38 * t17;
t49 = t1 * qJD(1);
t48 = qJD(1) * t30;
t47 = qJD(1) * t32;
t46 = qJD(2) * t30;
t45 = qJD(2) * t32;
t44 = t30 * qJD(4);
t27 = t32 * qJD(4);
t42 = t13 * t48;
t41 = t13 * t47;
t40 = t16 * t48;
t39 = pkin(1) * t43;
t37 = t56 - t61;
t2 = t37 * t30;
t35 = t2 * qJD(1) - t25 * t46;
t3 = t37 * t32;
t34 = t3 * qJD(1) - t25 * t45;
t24 = t28 * pkin(2) + pkin(6);
t23 = t30 * t27;
t21 = t22 * qJD(4);
t15 = t30 * t59;
t14 = pkin(6) + t38;
t10 = t16 * t46;
t5 = t60 * t32;
t4 = t60 * t30;
t6 = [0, 0, 0, 0, -t31 * t50, -t33 * t50, t1 * qJD(2), t23, t21, 0, 0, 0, t13 * t44 - t16 * t45, t13 * t27 + t10; 0, 0, 0, 0, -t31 * t39, -t33 * t39, t49 + (-t16 * t29 + t17 * t28) * qJD(2) * pkin(2), t23, t21, 0, 0, 0, t4 * qJD(4) - t16 * t59, t5 * qJD(4) + t10 + t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t15, t58, t27, -t44, 0, t4 * qJD(2) - t14 * t27 + t42, t5 * qJD(2) + t14 * t44 + t41; 0, 0, 0, 0, t31 * t51, t33 * t51, -t49, t23, t21, 0, 0, 0, -t2 * qJD(4) + t16 * t47, -t3 * qJD(4) - t40; 0, 0, 0, 0, 0, 0, 0, t23, t21, 0, 0, 0, t25 * t44, t25 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t15, t58, t27, -t44, 0, -t24 * t27 - t35, t24 * t44 - t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t27; 0, 0, 0, 0, 0, 0, 0, -t15, -t58, 0, 0, 0, t2 * qJD(2) - t42, t3 * qJD(2) - t41; 0, 0, 0, 0, 0, 0, 0, -t15, -t58, 0, 0, 0, t35, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t6;
