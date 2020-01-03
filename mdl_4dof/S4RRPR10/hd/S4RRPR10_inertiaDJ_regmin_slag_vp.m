% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x21]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRPR10_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR10_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_inertiaDJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:11:55
% EndTime: 2019-12-31 17:11:56
% DurationCPUTime: 0.28s
% Computational Cost: add. (138->57), mult. (396->126), div. (0->0), fcn. (256->4), ass. (0->55)
t24 = cos(qJ(2));
t20 = t24 ^ 2;
t22 = sin(qJ(2));
t36 = qJD(2) * (t22 ^ 2 - t20);
t21 = sin(qJ(4));
t17 = t21 ^ 2;
t23 = cos(qJ(4));
t56 = -t23 ^ 2 + t17;
t35 = t56 * qJD(4);
t25 = -pkin(2) - pkin(6);
t53 = t22 * qJ(3);
t29 = -t24 * t25 + t53;
t45 = t24 * qJD(3);
t57 = pkin(3) + pkin(5);
t12 = t57 * t24;
t47 = t12 * qJD(4);
t59 = t29 * qJD(2) - t45 - t47;
t58 = 0.2e1 * qJD(3);
t54 = qJ(3) * t24;
t52 = qJD(2) * t12;
t51 = qJD(4) * t21;
t50 = qJD(4) * t23;
t49 = qJD(4) * t24;
t48 = qJD(4) * t25;
t46 = t22 * qJD(2);
t16 = t24 * qJD(2);
t44 = qJ(3) * qJD(4);
t43 = -0.2e1 * pkin(1) * qJD(2);
t42 = pkin(5) * t46;
t41 = t21 * t49;
t40 = t23 * t49;
t39 = t22 * t16;
t38 = t23 * t46;
t37 = t21 * t50;
t34 = pkin(2) * t46 - t22 * qJD(3);
t33 = t21 * t38;
t11 = t57 * t22;
t7 = -pkin(1) - t29;
t32 = t23 * t11 - t21 * t7;
t31 = t21 * t11 + t23 * t7;
t30 = -t24 * pkin(2) - t53;
t8 = t57 * t46;
t27 = -t8 + (-t22 * t25 - t54) * qJD(4);
t26 = t30 * qJD(2) + t45;
t14 = pkin(5) * t16;
t13 = 0.2e1 * t39;
t10 = -pkin(1) + t30;
t9 = pkin(3) * t16 + t14;
t6 = -t21 * t16 - t22 * t50;
t5 = t23 * t16 - t22 * t51;
t4 = -qJ(3) * t16 + t34;
t3 = (pkin(6) * t22 - t54) * qJD(2) + t34;
t2 = -t31 * qJD(4) - t21 * t3 + t23 * t9;
t1 = -t32 * qJD(4) - t21 * t9 - t23 * t3;
t15 = [0, 0, 0, t13, -0.2e1 * t36, 0, 0, 0, t22 * t43, t24 * t43, 0, -0.2e1 * t10 * t46 + 0.2e1 * t4 * t24, -0.2e1 * t10 * t16 - 0.2e1 * t4 * t22, 0.2e1 * t10 * t4, -0.2e1 * t17 * t39 + 0.2e1 * t20 * t37, -0.2e1 * t20 * t35 - 0.4e1 * t24 * t33, 0.2e1 * t21 * t36 - 0.2e1 * t22 * t40, 0.2e1 * t22 * t41 + 0.2e1 * t23 * t36, t13, 0.2e1 * (-t23 * t52 + t2) * t22 + 0.2e1 * (t32 * qJD(2) - t21 * t47 - t8 * t23) * t24, 0.2e1 * (t21 * t52 + t1) * t22 + 0.2e1 * (-t31 * qJD(2) + t8 * t21 - t23 * t47) * t24; 0, 0, 0, 0, 0, t16, -t46, 0, -t14, t42, t26, t14, -t42, t26 * pkin(5), t24 * t35 + t33, 0.4e1 * t24 * t37 - t56 * t46, t5, t6, 0, t27 * t21 - t59 * t23, t59 * t21 + t27 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, qJ(3) * t58, -0.2e1 * t37, 0.2e1 * t35, 0, 0, 0, 0.2e1 * qJD(3) * t21 + 0.2e1 * t23 * t44, 0.2e1 * qJD(3) * t23 - 0.2e1 * t21 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, 0, t14, 0, 0, 0, 0, 0, t5, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21 * t46 - t40, t38 + t41, t16, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, -t50, 0, -t21 * t48, -t23 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, -t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t15;
