% Calculate inertial parameters regressor of joint inertia matrix for
% S5PRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRPR8_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t28 = sin(pkin(9));
t47 = t28 * pkin(3);
t22 = pkin(7) + t47;
t30 = sin(qJ(5));
t26 = t30 ^ 2;
t33 = cos(qJ(5));
t27 = t33 ^ 2;
t41 = t26 + t27;
t52 = t41 * t22;
t31 = sin(qJ(3));
t32 = sin(qJ(2));
t34 = cos(qJ(3));
t35 = cos(qJ(2));
t15 = -t31 * t32 + t34 * t35;
t16 = t31 * t35 + t34 * t32;
t29 = cos(pkin(9));
t4 = -t29 * t15 + t28 * t16;
t51 = t4 ^ 2;
t50 = 0.2e1 * t30;
t49 = -0.2e1 * t33;
t25 = t34 * pkin(2);
t24 = t25 + pkin(3);
t45 = t31 * pkin(2);
t40 = t29 * t45;
t13 = t28 * t24 + t40;
t11 = pkin(7) + t13;
t48 = t11 * t41;
t46 = t29 * pkin(3);
t44 = t4 * t33;
t38 = -t29 * t24 + t28 * t45;
t10 = -pkin(4) + t38;
t23 = -pkin(4) - t46;
t43 = t10 + t23;
t6 = t28 * t15 + t29 * t16;
t1 = t41 * t6;
t21 = t33 * t50;
t3 = t6 ^ 2;
t2 = t4 * t30;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32 ^ 2 + t35 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15 ^ 2 + t16 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3 + t51, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41 * t3 + t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t32, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t16, 0, (t15 * t34 + t16 * t31) * pkin(2), 0, 0, 0, 0, 0, 0, -t4, -t6, 0, t6 * t13 + t38 * t4, 0, 0, 0, 0, 0, 0, -t44, t2, t1, t11 * t1 + t4 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t25, -0.2e1 * t45, 0, (t31 ^ 2 + t34 ^ 2) * pkin(2) ^ 2, 0, 0, 0, 0, 0, 1, -0.2e1 * t38, -0.2e1 * t13, 0, t13 ^ 2 + t38 ^ 2, t26, t21, 0, t27, 0, 0, t10 * t49, t10 * t50, 0.2e1 * t48, t41 * t11 ^ 2 + t10 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t16, 0, 0, 0, 0, 0, 0, 0, 0, -t4, -t6, 0, (t28 * t6 - t29 * t4) * pkin(3), 0, 0, 0, 0, 0, 0, -t44, t2, t1, t4 * t23 + t52 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t25, -t45, 0, 0, 0, 0, 0, 0, 0, 1, -t38 + t46, -t40 + (-pkin(3) - t24) * t28, 0, (t13 * t28 - t29 * t38) * pkin(3), t26, t21, 0, t27, 0, 0, -t43 * t33, t43 * t30, t52 + t48, t10 * t23 + t11 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t46, -0.2e1 * t47, 0, (t28 ^ 2 + t29 ^ 2) * pkin(3) ^ 2, t26, t21, 0, t27, 0, 0, t23 * t49, t23 * t50, 0.2e1 * t52, t41 * t22 ^ 2 + t23 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30 * t6, -t33 * t6, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, t33, 0, -t30 * t11, -t33 * t11, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, t33, 0, -t30 * t22, -t33 * t22, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t30, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t5;
