% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:37
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPPR5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR5_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:36:30
% EndTime: 2021-01-15 19:36:32
% DurationCPUTime: 0.33s
% Computational Cost: add. (463->77), mult. (1095->143), div. (0->0), fcn. (945->6), ass. (0->53)
t65 = 2 * qJD(4);
t64 = -pkin(3) - pkin(4);
t63 = -qJ(3) - pkin(6);
t49 = sin(qJ(2));
t37 = t63 * t49;
t51 = cos(qJ(2));
t38 = t63 * t51;
t46 = sin(pkin(8));
t47 = cos(pkin(8));
t23 = t46 * t37 - t47 * t38;
t48 = sin(qJ(5));
t62 = qJD(5) * t48;
t50 = cos(qJ(5));
t61 = qJD(5) * t50;
t60 = t49 * qJD(2);
t59 = t51 * qJD(2);
t58 = -0.2e1 * pkin(1) * qJD(2);
t45 = pkin(2) * t60;
t44 = -t51 * pkin(2) - pkin(1);
t43 = -t47 * pkin(2) - pkin(3);
t56 = qJD(2) * t63;
t29 = t51 * qJD(3) + t49 * t56;
t53 = -t49 * qJD(3) + t51 * t56;
t13 = t46 * t29 - t47 * t53;
t14 = t47 * t29 + t46 * t53;
t22 = -t47 * t37 - t46 * t38;
t57 = t22 * t13 + t23 * t14;
t32 = t46 * t49 - t47 * t51;
t33 = t46 * t51 + t47 * t49;
t19 = t48 * t32 + t50 * t33;
t55 = t33 * qJ(4) - t44;
t31 = -t46 * t60 + t47 * t59;
t54 = t31 * qJ(4) + t33 * qJD(4) - t45;
t30 = t33 * qJD(2);
t52 = 0.2e1 * t13 * t33 - 0.2e1 * t14 * t32 + 0.2e1 * t22 * t31 - 0.2e1 * t23 * t30;
t41 = t46 * pkin(2) + qJ(4);
t40 = -pkin(4) + t43;
t21 = t48 * qJD(4) + (t40 * t48 + t41 * t50) * qJD(5);
t20 = -t50 * qJD(4) + (-t40 * t50 + t41 * t48) * qJD(5);
t18 = -t50 * t32 + t48 * t33;
t17 = t32 * pkin(3) - t55;
t16 = t32 * pkin(7) + t23;
t15 = -t33 * pkin(7) + t22;
t11 = t64 * t32 + t55;
t9 = t30 * pkin(3) - t54;
t8 = t30 * pkin(7) + t14;
t7 = -t31 * pkin(7) + t13;
t5 = t64 * t30 + t54;
t4 = t19 * qJD(5) - t50 * t30 + t48 * t31;
t3 = -t48 * t30 - t50 * t31 - t32 * t61 + t33 * t62;
t2 = t48 * t8 - t50 * t7 + (t15 * t48 + t16 * t50) * qJD(5);
t1 = -t48 * t7 - t50 * t8 + (-t15 * t50 + t16 * t48) * qJD(5);
t6 = [0, 0, 0, 0.2e1 * t49 * t59, 0.2e1 * (-t49 ^ 2 + t51 ^ 2) * qJD(2), 0, 0, 0, t49 * t58, t51 * t58, 0.2e1 * t44 * t30 + 0.2e1 * t32 * t45, 0.2e1 * t44 * t31 + 0.2e1 * t33 * t45, t52, 0.2e1 * t44 * t45 + 0.2e1 * t57, 0.2e1 * t17 * t30 + 0.2e1 * t9 * t32, t52, -0.2e1 * t17 * t31 - 0.2e1 * t9 * t33, 0.2e1 * t17 * t9 + 0.2e1 * t57, -0.2e1 * t19 * t3, 0.2e1 * t3 * t18 - 0.2e1 * t19 * t4, 0, 0, 0, 0.2e1 * t11 * t4 + 0.2e1 * t5 * t18, -0.2e1 * t11 * t3 + 0.2e1 * t5 * t19; 0, 0, 0, 0, 0, t59, -t60, 0, -pkin(6) * t59, pkin(6) * t60, -t13, -t14, (-t30 * t46 - t31 * t47) * pkin(2), (-t13 * t47 + t14 * t46) * pkin(2), -t13, -qJD(4) * t32 - t41 * t30 + t43 * t31, t14, t23 * qJD(4) + t13 * t43 + t14 * t41, 0, 0, t3, t4, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, t41 * t65, 0, 0, 0, 0, 0, 0.2e1 * t21, -0.2e1 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t31, 0, t45, t30, 0, -t31, t9, 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, t13, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, -t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t6;
