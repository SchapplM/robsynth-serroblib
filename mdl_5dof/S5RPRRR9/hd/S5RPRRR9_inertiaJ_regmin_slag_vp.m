% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRR9_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR9_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR9_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:07:59
% EndTime: 2019-12-31 19:08:01
% DurationCPUTime: 0.34s
% Computational Cost: add. (396->50), mult. (806->103), div. (0->0), fcn. (1010->8), ass. (0->51)
t31 = sin(pkin(9));
t32 = cos(pkin(9));
t35 = sin(qJ(3));
t37 = cos(qJ(3));
t19 = t35 * t31 - t37 * t32;
t24 = -t32 * pkin(2) - pkin(1);
t16 = t19 * pkin(3) + t24;
t56 = 0.2e1 * t16;
t55 = 0.2e1 * t24;
t34 = sin(qJ(4));
t54 = t34 * pkin(3);
t36 = cos(qJ(5));
t20 = t37 * t31 + t35 * t32;
t47 = pkin(6) + qJ(2);
t21 = t47 * t31;
t22 = t47 * t32;
t43 = -t37 * t21 - t35 * t22;
t39 = -t20 * pkin(7) + t43;
t51 = cos(qJ(4));
t40 = t35 * t21 - t37 * t22;
t9 = -t19 * pkin(7) - t40;
t4 = t34 * t9 - t51 * t39;
t53 = t4 * t36;
t44 = t51 * pkin(3);
t26 = -t44 - pkin(4);
t52 = pkin(4) - t26;
t14 = t51 * t19 + t34 * t20;
t33 = sin(qJ(5));
t11 = t33 * t14;
t15 = -t34 * t19 + t51 * t20;
t50 = t33 * t15;
t49 = t33 * t36;
t48 = t36 * t15;
t46 = t31 ^ 2 + t32 ^ 2;
t45 = -0.2e1 * t15 * t14;
t42 = -pkin(4) * t15 - pkin(8) * t14;
t25 = pkin(8) + t54;
t41 = -t14 * t25 + t15 * t26;
t30 = t36 ^ 2;
t29 = t33 ^ 2;
t23 = 0.2e1 * t49;
t13 = t15 ^ 2;
t12 = t36 * t14;
t10 = t33 * t48;
t7 = (-t29 + t30) * t15;
t6 = t14 * pkin(4) - t15 * pkin(8) + t16;
t5 = t34 * t39 + t51 * t9;
t3 = t4 * t33;
t2 = t33 * t6 + t36 * t5;
t1 = -t33 * t5 + t36 * t6;
t8 = [1, 0, 0, 0.2e1 * pkin(1) * t32, -0.2e1 * pkin(1) * t31, 0.2e1 * t46 * qJ(2), t46 * qJ(2) ^ 2 + pkin(1) ^ 2, t20 ^ 2, -0.2e1 * t20 * t19, 0, 0, 0, t19 * t55, t20 * t55, t13, t45, 0, 0, 0, t14 * t56, t15 * t56, t30 * t13, -0.2e1 * t13 * t49, 0.2e1 * t14 * t48, t33 * t45, t14 ^ 2, 0.2e1 * t1 * t14 + 0.2e1 * t4 * t50, -0.2e1 * t2 * t14 + 0.2e1 * t4 * t48; 0, 0, 0, -t32, t31, 0, -pkin(1), 0, 0, 0, 0, 0, t19, t20, 0, 0, 0, 0, 0, t14, t15, 0, 0, 0, 0, 0, t12, -t11; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, 0, t43, t40, 0, 0, t15, -t14, 0, -t4, -t5, t10, t7, t11, t12, 0, t41 * t33 - t53, t41 * t36 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t44, -0.2e1 * t54, t29, t23, 0, 0, 0, -0.2e1 * t26 * t36, 0.2e1 * t26 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t14, 0, -t4, -t5, t10, t7, t11, t12, 0, t42 * t33 - t53, t42 * t36 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t44, -t54, t29, t23, 0, 0, 0, t52 * t36, -t52 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t29, t23, 0, 0, 0, 0.2e1 * pkin(4) * t36, -0.2e1 * pkin(4) * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t50, t14, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, t36, 0, -t33 * t25, -t36 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, t36, 0, -t33 * pkin(8), -t36 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t8;
