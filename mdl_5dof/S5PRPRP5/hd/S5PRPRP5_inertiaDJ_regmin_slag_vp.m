% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x19]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRPRP5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP5_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:38:41
% EndTime: 2019-12-05 15:38:42
% DurationCPUTime: 0.30s
% Computational Cost: add. (257->61), mult. (732->116), div. (0->0), fcn. (665->6), ass. (0->41)
t27 = sin(pkin(8));
t28 = cos(pkin(8));
t51 = cos(qJ(4));
t39 = qJD(4) * t51;
t29 = sin(qJ(4));
t46 = qJD(4) * t29;
t54 = -t27 * t46 + t28 * t39;
t32 = -t29 * t27 + t51 * t28;
t53 = 0.2e1 * t54;
t52 = 2 * qJD(5);
t49 = t29 * t28;
t48 = pkin(6) + qJ(3);
t47 = t27 ^ 2 + t28 ^ 2;
t30 = sin(qJ(2));
t24 = t30 * qJD(2);
t31 = cos(qJ(2));
t45 = t31 * qJD(2);
t43 = t30 * t45;
t23 = -t28 * pkin(3) - pkin(2);
t41 = t48 * t27;
t40 = t47 * t31;
t38 = t51 * qJD(3);
t37 = t47 * qJD(3);
t36 = 0.2e1 * t37;
t34 = t51 * t41;
t16 = t51 * t27 + t49;
t33 = t16 * t24 - t31 * t54;
t13 = t16 * qJD(4);
t19 = t48 * t28;
t11 = t32 * t30;
t10 = t16 * t30;
t9 = t51 * t19 - t29 * t41;
t8 = t29 * t19 + t34;
t7 = -pkin(4) * t32 - t16 * qJ(5) + t23;
t6 = -t31 * t13 - t24 * t32;
t5 = t16 * t45 + t54 * t30;
t4 = t30 * t13 - t32 * t45;
t3 = t19 * t39 + qJD(3) * t49 + (-t48 * t46 + t38) * t27;
t2 = qJD(4) * t34 - t28 * t38 + (qJD(3) * t27 + qJD(4) * t19) * t29;
t1 = t13 * pkin(4) - qJ(5) * t54 - t16 * qJD(5);
t12 = [0, 0, 0, 0, 0, 0, 0, 0.2e1 * (-0.1e1 + t47) * t43, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t10 * t5 - 0.2e1 * t11 * t4 - 0.2e1 * t43; 0, 0, -t24, -t45, -t28 * t24, t27 * t24, qJD(2) * t40, t30 * t37 + (-pkin(2) * t30 + qJ(3) * t40) * qJD(2), 0, 0, 0, 0, 0, t6, t33, t6, t10 * t54 - t11 * t13 + t5 * t16 - t32 * t4, -t33, -t31 * t1 + t10 * t3 - t11 * t2 + t7 * t24 - t4 * t9 + t5 * t8; 0, 0, 0, 0, 0, 0, t36, qJ(3) * t36, t16 * t53, -0.2e1 * t16 * t13 + 0.2e1 * t32 * t54, 0, 0, 0, 0.2e1 * t23 * t13, t23 * t53, -0.2e1 * t1 * t32 + 0.2e1 * t7 * t13, -0.2e1 * t9 * t13 + 0.2e1 * t3 * t16 - 0.2e1 * t2 * t32 + 0.2e1 * t54 * t8, -0.2e1 * t1 * t16 - 0.2e1 * t54 * t7, 0.2e1 * t7 * t1 - 0.2e1 * t9 * t2 + 0.2e1 * t8 * t3; 0, 0, 0, 0, 0, 0, 0, t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t54, t13, 0, -t54, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, t4, -t5, 0, -t4, -t5 * pkin(4) - t4 * qJ(5) + t11 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, -t13, 0, -t3, t2, -t3, -pkin(4) * t54 - t13 * qJ(5) + qJD(5) * t32, -t2, -t3 * pkin(4) - t2 * qJ(5) + t9 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, qJ(5) * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t12;
