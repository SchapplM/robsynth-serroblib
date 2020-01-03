% Calculate minimal parameter regressor of coriolis matrix for
% S4RRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x18]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRPP5_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP5_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP5_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP5_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:00:38
% EndTime: 2019-12-31 17:00:39
% DurationCPUTime: 0.26s
% Computational Cost: add. (202->62), mult. (359->74), div. (0->0), fcn. (270->2), ass. (0->54)
t29 = cos(qJ(2));
t20 = t29 * qJD(3);
t28 = sin(qJ(2));
t52 = t28 * qJ(3);
t31 = -t29 * pkin(2) - t52;
t58 = t31 * qJD(2) + t20;
t27 = pkin(2) + qJ(4);
t57 = -t27 * t29 - t52;
t56 = pkin(3) + pkin(5);
t8 = -pkin(1) + t57;
t12 = t28 * pkin(2) - t29 * qJ(3);
t9 = t28 * qJ(4) + t12;
t1 = t8 * t9;
t54 = t1 * qJD(1);
t2 = t9 * t28 + t8 * t29;
t53 = t2 * qJD(1);
t3 = -t8 * t28 + t9 * t29;
t51 = t3 * qJD(1);
t11 = -pkin(1) + t31;
t4 = t11 * t29 + t12 * t28;
t50 = t4 * qJD(1);
t5 = -t11 * t28 + t12 * t29;
t49 = t5 * qJD(1);
t48 = qJD(1) * t28;
t47 = qJD(1) * t29;
t46 = qJD(2) * t28;
t21 = qJD(2) * t29;
t45 = qJD(3) * t28;
t44 = qJD(4) * t29;
t14 = t56 * t29;
t43 = t14 * qJD(2);
t25 = t28 ^ 2;
t26 = t29 ^ 2;
t15 = t26 - t25;
t42 = t15 * qJD(1);
t41 = t26 * qJD(1);
t40 = t27 * qJD(2);
t39 = pkin(1) * t48;
t38 = pkin(1) * t47;
t37 = pkin(5) * t46;
t36 = t8 * t48;
t35 = t8 * t47;
t34 = t11 * t12 * qJD(1);
t33 = t11 * t48;
t32 = t28 * t20;
t24 = qJ(3) * qJD(3);
t23 = qJD(2) * qJ(3);
t19 = t25 * qJD(1);
t18 = t25 * qJD(3);
t17 = pkin(5) * t21;
t16 = t28 * t47;
t13 = t56 * t28;
t10 = t13 * qJD(2);
t6 = [0, 0, 0, t28 * t21, t15 * qJD(2), 0, 0, 0, -pkin(1) * t46, -pkin(1) * t21, 0, t5 * qJD(2) - t32, -t4 * qJD(2) + t18, (qJD(2) * t12 - t45) * t11, 0, -t2 * qJD(2) + t28 * t44 + t18, -t3 * qJD(2) + t26 * qJD(4) + t32, t1 * qJD(2) + (-t44 - t45) * t8; 0, 0, 0, t16, t42, t21, -t46, 0, -t17 - t39, t37 - t38, t58, t17 + t49, -t37 - t50, t58 * pkin(5) + t34, t57 * qJD(2) - qJD(4) * t28 + t20, -t10 - t53, -t43 - t51, t54 + (-t13 * qJ(3) - t14 * t27) * qJD(2) + t14 * qJD(3) - t13 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t16, t19, t17 - t33, t21, t19, t16, -t36 + t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, t16, t41, -t10 - t35; 0, 0, 0, -t16, -t42, 0, 0, 0, t39, t38, 0, -t49, t50, -t34, 0, t53, t51, -t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t24, 0, qJD(3), qJD(4), t27 * qJD(4) + t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), t23, 0, qJD(2), 0, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t19, t33, 0, -t19, -t16, t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2), -t23, 0, -qJD(2), 0, -qJD(4) - t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t41, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2), qJD(3) - t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t6;
