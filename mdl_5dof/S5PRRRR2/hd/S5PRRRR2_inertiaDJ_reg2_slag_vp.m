% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRRR2_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR2_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_inertiaDJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:04:51
% EndTime: 2019-12-05 17:04:54
% DurationCPUTime: 0.53s
% Computational Cost: add. (234->47), mult. (774->87), div. (0->0), fcn. (480->6), ass. (0->48)
t32 = sin(qJ(5));
t30 = t32 ^ 2;
t35 = cos(qJ(5));
t31 = t35 ^ 2;
t53 = t30 + t31;
t36 = cos(qJ(4));
t56 = cos(qJ(3));
t45 = t56 * pkin(2);
t39 = t45 + pkin(3);
t38 = t36 * t39;
t40 = qJD(3) * t45;
t33 = sin(qJ(4));
t34 = sin(qJ(3));
t60 = pkin(2) * t34;
t49 = t33 * t60;
t61 = qJD(3) + qJD(4);
t7 = -qJD(4) * t38 - t36 * t40 + t61 * t49;
t43 = t53 * t7;
t54 = t34 * t36;
t63 = t61 * pkin(2) * (t56 * t33 + t54);
t15 = -t38 + t49;
t51 = qJD(4) * t33;
t47 = pkin(3) * t51;
t8 = t47 + t63;
t59 = t15 * t8;
t58 = t33 * pkin(3);
t57 = t36 * t8;
t29 = t35 * qJD(5);
t2 = t15 * t29 + t8 * t32;
t55 = t15 * t33;
t52 = qJD(4) * pkin(3);
t46 = t36 * t52;
t12 = t53 * t46;
t16 = pkin(2) * t54 + t33 * t39;
t50 = t32 * qJD(5);
t48 = qJD(3) * t60;
t44 = t32 * t29;
t3 = t15 * t50 - t8 * t35;
t42 = t53 * t36;
t41 = -0.2e1 * t47;
t13 = -pkin(3) * t36 * t29 + t32 * t47;
t14 = (-t35 * t51 - t36 * t50) * pkin(3);
t28 = pkin(6) + t58;
t26 = -0.2e1 * t44;
t25 = 0.2e1 * t44;
t17 = 0.2e1 * (-t30 + t31) * qJD(5);
t11 = pkin(6) + t16;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t48, -0.2e1 * t40, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t8, 0.2e1 * t7, 0, -0.2e1 * t16 * t7 + 0.2e1 * t59, t25, t17, 0, t26, 0, 0, 0.2e1 * t3, 0.2e1 * t2, -0.2e1 * t43, -0.2e1 * t11 * t43 + 0.2e1 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, -t40, 0, 0, 0, 0, 0, 0, 0, 0, t41 - t63, t7 - t46, 0, (-t33 * t7 - t57 + (t16 * t36 + t55) * qJD(4)) * pkin(3), t25, t17, 0, t26, 0, 0, t14 + t3, t13 + t2, t12 - t43, -t28 * t43 + (-t57 + (t11 * t42 + t55) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -0.2e1 * t46, 0, 0, t25, t17, 0, t26, 0, 0, 0.2e1 * t14, 0.2e1 * t13, 0.2e1 * t12, 0.2e1 * (t53 * t28 - t58) * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, t7, 0, 0, t25, t17, 0, t26, 0, 0, t3, t2, -t43, -pkin(6) * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, -t46, 0, 0, t25, t17, 0, t26, 0, 0, t14, t13, t12, pkin(6) * t42 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t17, 0, t26, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, -t29, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, -t50, 0, -t11 * t29 + t32 * t7, t11 * t50 + t35 * t7, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, -t50, 0, -t28 * t29 - t32 * t46, t28 * t50 - t35 * t46, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, -t50, 0, -pkin(6) * t29, pkin(6) * t50, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
