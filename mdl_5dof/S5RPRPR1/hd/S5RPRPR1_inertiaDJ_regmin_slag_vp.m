% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x24]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:34
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR1_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:33:52
% EndTime: 2021-01-15 11:33:55
% DurationCPUTime: 0.34s
% Computational Cost: add. (450->64), mult. (943->124), div. (0->0), fcn. (849->6), ass. (0->46)
t44 = sin(pkin(8));
t45 = cos(pkin(8));
t49 = cos(qJ(3));
t57 = t49 * qJD(3);
t47 = sin(qJ(3));
t58 = t47 * qJD(3);
t26 = t44 * t58 - t45 * t57;
t27 = -t44 * t57 - t45 * t58;
t62 = (t26 * t44 - t27 * t45) * pkin(3);
t31 = -t44 * t47 + t45 * t49;
t32 = t44 * t49 + t45 * t47;
t46 = sin(qJ(5));
t48 = cos(qJ(5));
t14 = t48 * t31 - t46 * t32;
t5 = t14 * qJD(5) - t48 * t26 + t46 * t27;
t61 = 2 * qJD(2);
t60 = pkin(3) * t44;
t50 = -pkin(1) - pkin(6);
t59 = qJ(4) - t50;
t33 = t59 * t47;
t34 = t59 * t49;
t16 = -t45 * t33 - t44 * t34;
t40 = t47 * pkin(3) + qJ(2);
t35 = pkin(3) * t57 + qJD(2);
t56 = qJ(2) * qJD(3);
t24 = -t49 * qJD(4) + t59 * t58;
t25 = -qJD(3) * t34 - t47 * qJD(4);
t9 = t45 * t24 - t44 * t25;
t15 = t44 * t33 - t45 * t34;
t10 = t44 * t24 + t45 * t25;
t55 = t32 * t26 - t31 * t27;
t13 = t46 * t31 + t48 * t32;
t52 = t10 * t32 + t15 * t27 - t16 * t26 + t9 * t31;
t51 = -t13 * qJD(5) + t46 * t26 + t48 * t27;
t41 = t45 * pkin(3) + pkin(4);
t21 = (-t41 * t46 - t48 * t60) * qJD(5);
t20 = (-t41 * t48 + t46 * t60) * qJD(5);
t19 = t32 * pkin(4) + t40;
t17 = -t26 * pkin(4) + t35;
t12 = -t32 * pkin(7) + t16;
t11 = -t31 * pkin(7) + t15;
t8 = t26 * pkin(7) + t10;
t7 = -t27 * pkin(7) + t9;
t2 = -t46 * t8 + t48 * t7 + (-t11 * t46 - t12 * t48) * qJD(5);
t1 = -t46 * t7 - t48 * t8 + (-t11 * t48 + t12 * t46) * qJD(5);
t3 = [0, 0, 0, 0, t61, qJ(2) * t61, -0.2e1 * t47 * t57, 0.2e1 * (t47 ^ 2 - t49 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * qJD(2) * t47 + 0.2e1 * t49 * t56, 0.2e1 * qJD(2) * t49 - 0.2e1 * t47 * t56, -0.2e1 * t40 * t26 + 0.2e1 * t35 * t32, 0.2e1 * t40 * t27 + 0.2e1 * t35 * t31, -0.2e1 * t52, 0.2e1 * t16 * t10 + 0.2e1 * t15 * t9 + 0.2e1 * t40 * t35, 0.2e1 * t14 * t51, -0.2e1 * t13 * t51 - 0.2e1 * t14 * t5, 0, 0, 0, 0.2e1 * t17 * t13 + 0.2e1 * t19 * t5, 0.2e1 * t17 * t14 + 0.2e1 * t19 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t55, t52, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t55, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t58, -t57, 0, -t50 * t58, -t50 * t57, t9, -t10, t62, (t10 * t44 + t45 * t9) * pkin(3), 0, 0, t51, -t5, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, -t57, t27, t26, 0, -t62, 0, 0, 0, 0, 0, t51, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t21, 0.2e1 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, t27, 0, t35, 0, 0, 0, 0, 0, t5, t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, -t5, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t3;
