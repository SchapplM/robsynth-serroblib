% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:12
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPR5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 23:10:44
% EndTime: 2021-01-15 23:10:46
% DurationCPUTime: 0.42s
% Computational Cost: add. (556->65), mult. (1098->135), div. (0->0), fcn. (1299->8), ass. (0->60)
t39 = sin(qJ(3));
t40 = sin(qJ(2));
t42 = cos(qJ(3));
t43 = cos(qJ(2));
t25 = t39 * t40 - t42 * t43;
t32 = -pkin(2) * t43 - pkin(1);
t19 = pkin(3) * t25 + t32;
t65 = 0.2e1 * t19;
t64 = 0.2e1 * t32;
t38 = sin(qJ(5));
t63 = 0.2e1 * t38;
t41 = cos(qJ(5));
t62 = -0.2e1 * t41;
t61 = 0.2e1 * t43;
t60 = pkin(6) + pkin(7);
t36 = sin(pkin(9));
t59 = t36 * pkin(3);
t37 = cos(pkin(9));
t58 = t37 * pkin(3);
t57 = t39 * pkin(2);
t49 = t60 * t43;
t50 = t60 * t40;
t18 = t39 * t50 - t42 * t49;
t10 = -t25 * qJ(4) - t18;
t17 = -t39 * t49 - t42 * t50;
t26 = t39 * t43 + t40 * t42;
t45 = -t26 * qJ(4) + t17;
t5 = t36 * t10 - t37 * t45;
t56 = t5 * t41;
t15 = t37 * t25 + t26 * t36;
t12 = t38 * t15;
t16 = -t25 * t36 + t26 * t37;
t55 = t38 * t16;
t54 = t38 * t41;
t53 = t41 * t16;
t33 = t42 * pkin(2);
t31 = t33 + pkin(3);
t48 = -t37 * t31 + t36 * t57;
t20 = -pkin(4) + t48;
t30 = -pkin(4) - t58;
t52 = t20 + t30;
t51 = t37 * t57;
t23 = t31 * t36 + t51;
t21 = pkin(8) + t23;
t47 = -t15 * t21 + t16 * t20;
t29 = pkin(8) + t59;
t46 = -t15 * t29 + t16 * t30;
t35 = t41 ^ 2;
t34 = t38 ^ 2;
t28 = 0.2e1 * t54;
t14 = t16 ^ 2;
t13 = t41 * t15;
t11 = t38 * t53;
t8 = (-t34 + t35) * t16;
t7 = t37 * t10 + t36 * t45;
t4 = pkin(4) * t15 - pkin(8) * t16 + t19;
t3 = t5 * t38;
t2 = t38 * t4 + t41 * t7;
t1 = -t38 * t7 + t4 * t41;
t6 = [1, 0, 0, t40 ^ 2, t40 * t61, 0, 0, 0, pkin(1) * t61, -0.2e1 * pkin(1) * t40, t26 ^ 2, -0.2e1 * t26 * t25, 0, 0, 0, t25 * t64, t26 * t64, t15 * t65, t16 * t65, -0.2e1 * t15 * t7 + 0.2e1 * t16 * t5, t19 ^ 2 + t5 ^ 2 + t7 ^ 2, t35 * t14, -0.2e1 * t14 * t54, 0.2e1 * t15 * t53, -0.2e1 * t15 * t55, t15 ^ 2, 0.2e1 * t1 * t15 + 0.2e1 * t5 * t55, -0.2e1 * t15 * t2 + 0.2e1 * t5 * t53; 0, 0, 0, 0, 0, t40, t43, 0, -t40 * pkin(6), -t43 * pkin(6), 0, 0, t26, -t25, 0, t17, t18, -t5, -t7, -t15 * t23 + t16 * t48, t23 * t7 + t48 * t5, t11, t8, t12, t13, 0, t38 * t47 - t56, t41 * t47 + t3; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t33, -0.2e1 * t57, -0.2e1 * t48, -0.2e1 * t23, 0, t23 ^ 2 + t48 ^ 2, t34, t28, 0, 0, 0, t20 * t62, t20 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t25, 0, t17, t18, -t5, -t7, (-t15 * t36 - t16 * t37) * pkin(3), (t36 * t7 - t37 * t5) * pkin(3), t11, t8, t12, t13, 0, t38 * t46 - t56, t41 * t46 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t33, -t57, -t48 + t58, -t51 + (-pkin(3) - t31) * t36, 0, (t23 * t36 - t37 * t48) * pkin(3), t34, t28, 0, 0, 0, -t52 * t41, t52 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * t58, -0.2e1 * t59, 0, (t36 ^ 2 + t37 ^ 2) * pkin(3) ^ 2, t34, t28, 0, 0, 0, t30 * t62, t30 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, t16, 0, t19, 0, 0, 0, 0, 0, t13, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, -t55, t15, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t41, 0, -t38 * t21, -t41 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t41, 0, -t38 * t29, -t41 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t6;
