% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRPR9_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:24:37
% EndTime: 2019-12-31 18:24:40
% DurationCPUTime: 0.56s
% Computational Cost: add. (225->55), mult. (396->94), div. (0->0), fcn. (357->6), ass. (0->52)
t33 = sin(qJ(3));
t26 = t33 ^ 2;
t35 = cos(qJ(3));
t28 = t35 ^ 2;
t16 = t26 + t28;
t36 = -pkin(3) - pkin(7);
t56 = t36 * t35;
t55 = 0.2e1 * t35;
t54 = 2 * qJ(4);
t30 = sin(pkin(8));
t53 = t30 * pkin(1);
t31 = cos(pkin(8));
t52 = t31 * pkin(1);
t51 = t35 * pkin(3);
t32 = sin(qJ(5));
t50 = t32 * t33;
t49 = t32 * t35;
t48 = t33 * t35;
t34 = cos(qJ(5));
t47 = t34 * t32;
t46 = t34 * t35;
t19 = pkin(6) + t53;
t45 = t16 * t19 ^ 2;
t25 = t32 ^ 2;
t27 = t34 ^ 2;
t15 = t25 + t27;
t44 = t35 * qJ(4);
t43 = -0.2e1 * t48;
t42 = t32 * t46;
t20 = -pkin(2) - t52;
t24 = t33 * qJ(4);
t41 = t20 - t24;
t5 = t41 + t56;
t12 = t33 * t19;
t7 = t33 * pkin(4) + t12;
t2 = -t32 * t5 + t34 * t7;
t3 = t32 * t7 + t34 * t5;
t1 = t2 * t34 + t3 * t32;
t40 = -t33 * pkin(3) + t44;
t39 = t33 * t36 + t44;
t37 = qJ(4) ^ 2;
t23 = t34 * t33;
t22 = t27 * t28;
t21 = t25 * t28;
t18 = 0.2e1 * t48;
t14 = t35 * t19;
t10 = t15 * t36;
t9 = t15 * t35;
t8 = t35 * pkin(4) + t14;
t6 = t41 - t51;
t4 = 0.2e1 * t16 * t19;
t11 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t52, -0.2e1 * t53, 0, (t30 ^ 2 + t31 ^ 2) * pkin(1) ^ 2, t26, t18, 0, t28, 0, 0, -0.2e1 * t20 * t35, 0.2e1 * t20 * t33, t4, t20 ^ 2 + t45, 0, 0, 0, t26, t18, t28, t4, t6 * t55, -0.2e1 * t6 * t33, t6 ^ 2 + t45, t21, 0.2e1 * t28 * t47, t32 * t43, t22, t34 * t43, t26, 0.2e1 * t2 * t33 + 0.2e1 * t8 * t46, -0.2e1 * t3 * t33 - 0.2e1 * t8 * t49, (t2 * t32 - t3 * t34) * t55, t2 ^ 2 + t3 ^ 2 + t8 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * t35 + t8 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21 + t22 + t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, t35, 0, -t12, -t14, 0, 0, 0, -t33, -t35, 0, 0, 0, t40, t12, t14, t40 * t19, -t42, (t25 - t27) * t35, t23, t42, -t50, 0, t8 * t32 + t39 * t34, -t39 * t32 + t8 * t34, -t1, t8 * qJ(4) + t1 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t33, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, t33, t24 + t51, 0, 0, 0, 0, 0, 0, t50, t23, t9, -t15 * t56 + t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -0.2e1 * pkin(3), t54, pkin(3) ^ 2 + t37, t27, -0.2e1 * t47, 0, t25, 0, 0, t32 * t54, t34 * t54, -0.2e1 * t10, t15 * t36 ^ 2 + t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, 0, t12, 0, 0, 0, 0, 0, 0, t23, -t50, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(3), 0, 0, 0, 0, 0, 0, 0, 0, -t15, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, 0, -t46, t33, t2, -t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, t49, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, -t32, 0, t34 * t36, -t32 * t36, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t32, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t11;
