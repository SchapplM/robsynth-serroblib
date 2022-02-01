% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:20
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRP1_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:20:00
% EndTime: 2022-01-20 10:20:02
% DurationCPUTime: 0.57s
% Computational Cost: add. (378->86), mult. (960->141), div. (0->0), fcn. (656->6), ass. (0->58)
t46 = sin(pkin(8));
t51 = cos(qJ(2));
t62 = pkin(1) * qJD(2);
t47 = cos(pkin(8));
t49 = sin(qJ(2));
t66 = t47 * t49;
t21 = (t46 * t51 + t66) * t62;
t48 = sin(qJ(4));
t40 = t48 * qJD(4);
t39 = pkin(4) * t40;
t11 = t39 + t21;
t37 = t51 * pkin(1) + pkin(2);
t52 = -t46 * t49 * pkin(1) + t47 * t37;
t19 = -pkin(3) - t52;
t50 = cos(qJ(4));
t67 = t50 * pkin(4);
t14 = t19 - t67;
t42 = t50 * qJD(4);
t68 = t11 * t48 + t14 * t42;
t65 = t19 * t42 + t21 * t48;
t64 = pkin(1) * t66 + t46 * t37;
t36 = -t47 * pkin(2) - pkin(3);
t28 = t36 - t67;
t44 = t48 ^ 2;
t63 = qJD(4) * t44 * pkin(4) + t28 * t42;
t20 = pkin(7) + t64;
t61 = qJ(5) + t20;
t35 = t46 * pkin(2) + pkin(7);
t60 = qJ(5) + t35;
t59 = t49 * t62;
t58 = t51 * t62;
t57 = pkin(4) * t42;
t5 = t14 * t40;
t23 = t28 * t40;
t56 = t36 * t40;
t55 = t36 * t42;
t54 = t48 * t42;
t22 = -t46 * t59 + t47 * t58;
t45 = t50 ^ 2;
t4 = (t44 + t45) * t22;
t53 = t19 * t40 - t21 * t50;
t43 = t50 * qJ(5);
t41 = t50 * qJD(5);
t32 = -0.2e1 * t54;
t31 = 0.2e1 * t54;
t27 = 0.2e1 * (-t44 + t45) * qJD(4);
t26 = t50 * t35 + t43;
t25 = t60 * t48;
t17 = t50 * t22;
t16 = -t48 * qJD(5) - t60 * t42;
t15 = t60 * t40 - t41;
t10 = t15 * t50;
t8 = t50 * t20 + t43;
t7 = t61 * t48;
t3 = (-qJD(5) - t22) * t48 - t61 * t42;
t2 = t61 * t40 - t17 - t41;
t1 = t2 * t50;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t59, -0.2e1 * t58, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t21, -0.2e1 * t22, 0, -0.2e1 * t52 * t21 + 0.2e1 * t64 * t22, t31, t27, 0, t32, 0, 0, 0.2e1 * t53, 0.2e1 * t65, 0.2e1 * t4, 0.2e1 * t19 * t21 + 0.2e1 * t20 * t4, t31, t27, 0, t32, 0, 0, -0.2e1 * t11 * t50 + 0.2e1 * t5, 0.2e1 * t68, -0.2e1 * t3 * t48 - 0.2e1 * t1 + 0.2e1 * (-t48 * t8 + t50 * t7) * qJD(4), 0.2e1 * t14 * t11 - 0.2e1 * t8 * t2 - 0.2e1 * t7 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, -t58, 0, 0, 0, 0, 0, 0, 0, 0, -t21, -t22, 0, (-t21 * t47 + t22 * t46) * pkin(2), t31, t27, 0, t32, 0, 0, t53 + t56, t55 + t65, t4, t21 * t36 + t35 * t4, t31, t27, 0, t32, 0, 0, t23 + t5 + (-t11 - t39) * t50, t63 + t68, -t1 - t10 + (-t16 - t3) * t48 + ((t25 + t7) * t50 + (-t26 - t8) * t48) * qJD(4), pkin(4) * t5 + t11 * t28 - t8 * t15 - t7 * t16 - t2 * t26 - t3 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t27, 0, t32, 0, 0, 0.2e1 * t56, 0.2e1 * t55, 0, 0, t31, t27, 0, t32, 0, 0, -0.2e1 * pkin(4) * t54 + 0.2e1 * t23, 0.2e1 * t63, -0.2e1 * t16 * t48 - 0.2e1 * t10 + 0.2e1 * (t25 * t50 - t26 * t48) * qJD(4), 0.2e1 * pkin(4) * t23 - 0.2e1 * t26 * t15 - 0.2e1 * t25 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2 * t48 + t3 * t50 + (t48 * t7 + t50 * t8) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15 * t48 + t16 * t50 + (t25 * t48 + t26 * t50) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, -t40, 0, -t20 * t42 - t48 * t22, t20 * t40 - t17, 0, 0, 0, 0, t42, 0, -t40, 0, t3, t2, -t57, t3 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, -t40, 0, -t35 * t42, t35 * t40, 0, 0, 0, 0, t42, 0, -t40, 0, t16, t15, -t57, t16 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, -t42, 0, 0, 0, 0, 0, 0, 0, 0, -t40, -t42, 0, -t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t42, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t42, 0, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t6;
