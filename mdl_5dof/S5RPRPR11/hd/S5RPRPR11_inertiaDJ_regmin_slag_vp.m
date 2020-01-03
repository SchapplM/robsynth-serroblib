% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPR11_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR11_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:27:55
% EndTime: 2019-12-31 18:27:56
% DurationCPUTime: 0.30s
% Computational Cost: add. (392->70), mult. (940->124), div. (0->0), fcn. (877->6), ass. (0->46)
t41 = sin(pkin(8));
t56 = pkin(6) + qJ(2);
t33 = t56 * t41;
t42 = cos(pkin(8));
t34 = t56 * t42;
t44 = sin(qJ(3));
t58 = cos(qJ(3));
t18 = -t44 * t33 + t58 * t34;
t51 = qJD(3) * t58;
t52 = t58 * t42;
t10 = (qJD(2) * t41 + qJD(3) * t34) * t44 - qJD(2) * t52 + t33 * t51;
t57 = t44 * t41;
t25 = qJD(3) * t57 - t42 * t51;
t60 = -0.2e1 * t25;
t59 = 2 * qJD(4);
t46 = -pkin(3) - pkin(4);
t43 = sin(qJ(5));
t54 = qJD(5) * t43;
t45 = cos(qJ(5));
t53 = qJD(5) * t45;
t38 = -t42 * pkin(2) - pkin(1);
t50 = 0.2e1 * (t41 ^ 2 + t42 ^ 2) * qJD(2);
t17 = t58 * t33 + t44 * t34;
t31 = -t52 + t57;
t32 = t58 * t41 + t44 * t42;
t16 = t43 * t31 + t45 * t32;
t49 = -t25 * qJ(4) + t32 * qJD(4);
t47 = t32 * qJ(4) - t38;
t11 = t32 * qJD(2) + t18 * qJD(3);
t26 = t32 * qJD(3);
t23 = t43 * qJD(4) + (qJ(4) * t45 + t43 * t46) * qJD(5);
t22 = -t45 * qJD(4) + (qJ(4) * t43 - t45 * t46) * qJD(5);
t15 = -t45 * t31 + t43 * t32;
t14 = t31 * pkin(3) - t47;
t13 = t31 * pkin(7) + t18;
t12 = -t32 * pkin(7) + t17;
t9 = t46 * t31 + t47;
t8 = t26 * pkin(3) - t49;
t7 = t25 * pkin(7) + t11;
t6 = t26 * pkin(7) - t10;
t5 = t46 * t26 + t49;
t4 = t16 * qJD(5) - t43 * t25 - t45 * t26;
t3 = t45 * t25 - t43 * t26 - t31 * t53 + t32 * t54;
t2 = t43 * t6 - t45 * t7 + (t12 * t43 + t13 * t45) * qJD(5);
t1 = -t43 * t7 - t45 * t6 + (-t12 * t45 + t13 * t43) * qJD(5);
t19 = [0, 0, 0, 0, 0, t50, qJ(2) * t50, t32 * t60, 0.2e1 * t25 * t31 - 0.2e1 * t32 * t26, 0, 0, 0, 0.2e1 * t38 * t26, t38 * t60, 0.2e1 * t14 * t26 + 0.2e1 * t8 * t31, 0.2e1 * t10 * t31 + 0.2e1 * t11 * t32 - 0.2e1 * t17 * t25 - 0.2e1 * t18 * t26, 0.2e1 * t14 * t25 - 0.2e1 * t8 * t32, -0.2e1 * t18 * t10 + 0.2e1 * t17 * t11 + 0.2e1 * t14 * t8, -0.2e1 * t16 * t3, 0.2e1 * t3 * t15 - 0.2e1 * t16 * t4, 0, 0, 0, 0.2e1 * t5 * t15 + 0.2e1 * t9 * t4, 0.2e1 * t5 * t16 - 0.2e1 * t9 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t25, t26, 0, t25, t8, 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, -t26, 0, -t11, t10, -t11, pkin(3) * t25 - t26 * qJ(4) - t31 * qJD(4), -t10, -t11 * pkin(3) - t10 * qJ(4) + t18 * qJD(4), 0, 0, t3, t4, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, qJ(4) * t59, 0, 0, 0, 0, 0, 0.2e1 * t23, -0.2e1 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, 0, t11, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, -t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t19;
