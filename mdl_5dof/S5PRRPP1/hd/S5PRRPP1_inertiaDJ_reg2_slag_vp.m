% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRPP1_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP1_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:06:53
% EndTime: 2019-12-05 16:06:55
% DurationCPUTime: 0.44s
% Computational Cost: add. (325->55), mult. (847->109), div. (0->0), fcn. (711->4), ass. (0->42)
t32 = sin(pkin(8));
t34 = cos(qJ(3));
t33 = sin(qJ(3));
t49 = cos(pkin(8));
t41 = t49 * t33;
t23 = t32 * t34 + t41;
t20 = t23 * qJD(3);
t40 = t49 * t34;
t47 = t33 * qJD(3);
t21 = qJD(3) * t40 - t32 * t47;
t53 = t32 * t33;
t22 = -t40 + t53;
t51 = t20 * t23 + t22 * t21;
t57 = 2 * qJD(5);
t55 = t22 * t20;
t52 = -qJ(4) - pkin(6);
t48 = t23 * qJD(5);
t46 = t34 * qJD(3);
t45 = 0.2e1 * t55;
t44 = -0.2e1 * pkin(2) * qJD(3);
t31 = pkin(3) * t47;
t43 = t33 * t46;
t30 = -t34 * pkin(3) - pkin(2);
t25 = t52 * t34;
t16 = -t32 * t25 - t52 * t41;
t17 = -t49 * t25 + t52 * t53;
t39 = qJD(3) * t52;
t19 = t34 * qJD(4) + t33 * t39;
t36 = -t33 * qJD(4) + t34 * t39;
t7 = t32 * t19 - t49 * t36;
t8 = t49 * t19 + t32 * t36;
t42 = t16 * t7 + t17 * t8;
t15 = t23 * t21;
t38 = 0.2e1 * t15 + 0.2e1 * t55;
t37 = t20 * t16 + t21 * t17 + t22 * t7 + t23 * t8;
t35 = 0.2e1 * t16 * t21 - 0.2e1 * t17 * t20 - 0.2e1 * t8 * t22 + 0.2e1 * t7 * t23;
t29 = -t49 * pkin(3) - pkin(4);
t27 = t32 * pkin(3) + qJ(5);
t10 = 0.2e1 * t15;
t9 = t22 * pkin(4) - t23 * qJ(5) + t30;
t2 = t20 * pkin(4) - t21 * qJ(5) + t31 - t48;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t43, 0.2e1 * (-t33 ^ 2 + t34 ^ 2) * qJD(3), 0, -0.2e1 * t43, 0, 0, t33 * t44, t34 * t44, 0, 0, t10, -0.2e1 * t51, 0, t45, 0, 0, 0.2e1 * t30 * t20 + 0.2e1 * t22 * t31, 0.2e1 * t30 * t21 + 0.2e1 * t23 * t31, t35, 0.2e1 * t30 * t31 + 0.2e1 * t42, t10, 0, 0.2e1 * t51, 0, 0, t45, 0.2e1 * t2 * t22 + 0.2e1 * t9 * t20, t35, -0.2e1 * t2 * t23 - 0.2e1 * t9 * t21, 0.2e1 * t9 * t2 + 0.2e1 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, -t46, 0, 0, 0, 0, 0, 0, 0, 0, -t20, -t21, 0, (-t20 * t49 + t21 * t32) * pkin(3), 0, 0, 0, 0, 0, 0, -t20, 0, t21, t20 * t29 + t21 * t27 + t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0, -t47, 0, -pkin(6) * t46, pkin(6) * t47, 0, 0, 0, 0, t21, 0, -t20, 0, -t7, -t8, (-t20 * t32 - t21 * t49) * pkin(3), (t32 * t8 - t49 * t7) * pkin(3), 0, t21, 0, 0, t20, 0, -t7, -qJD(5) * t22 - t27 * t20 + t29 * t21, t8, t17 * qJD(5) + t8 * t27 + t7 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, t27 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, t21, 0, t31, 0, 0, 0, 0, 0, 0, t20, 0, -t21, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
