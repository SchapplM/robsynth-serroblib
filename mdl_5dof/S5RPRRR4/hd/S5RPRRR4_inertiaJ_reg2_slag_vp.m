% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPRRR4
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
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRR4_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t30 = sin(qJ(5));
t26 = t30 ^ 2;
t33 = cos(qJ(5));
t27 = t33 ^ 2;
t41 = t26 + t27;
t29 = cos(pkin(9));
t48 = t29 * pkin(1);
t19 = pkin(2) + t48;
t32 = sin(qJ(3));
t35 = cos(qJ(3));
t28 = sin(pkin(9));
t49 = t28 * pkin(1);
t11 = t35 * t19 - t32 * t49;
t10 = pkin(3) + t11;
t31 = sin(qJ(4));
t12 = t32 * t19 + t35 * t49;
t34 = cos(qJ(4));
t44 = t34 * t12;
t8 = t31 * t10 + t44;
t6 = pkin(8) + t8;
t51 = t41 * t6;
t47 = t31 * pkin(3);
t20 = pkin(8) + t47;
t52 = t41 * t20;
t50 = pkin(4) * t30;
t40 = -t34 * t10 + t31 * t12;
t5 = -pkin(4) + t40;
t46 = t5 * t33;
t24 = t34 * pkin(3);
t21 = -t24 - pkin(4);
t45 = t21 * t33;
t42 = pkin(8) * t41;
t25 = pkin(4) * t33;
t17 = 0.2e1 * t30 * t33;
t16 = t21 * t30;
t3 = t5 * t30;
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t48, -0.2e1 * t49, 0, (t28 ^ 2 + t29 ^ 2) * pkin(1) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t11, -0.2e1 * t12, 0, t11 ^ 2 + t12 ^ 2, 0, 0, 0, 0, 0, 1, -0.2e1 * t40, -0.2e1 * t8, 0, t40 ^ 2 + t8 ^ 2, t26, t17, 0, t27, 0, 0, -0.2e1 * t46, 0.2e1 * t3, 0.2e1 * t51, t41 * t6 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t11, -t12, 0, 0, 0, 0, 0, 0, 0, 1, t24 - t40, -t44 + (-pkin(3) - t10) * t31, 0, (t31 * t8 - t34 * t40) * pkin(3), t26, t17, 0, t27, 0, 0, (-t21 - t5) * t33, t16 + t3, t52 + t51, t5 * t21 + t52 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t24, -0.2e1 * t47, 0, (t31 ^ 2 + t34 ^ 2) * pkin(3) ^ 2, t26, t17, 0, t27, 0, 0, -0.2e1 * t45, 0.2e1 * t16, 0.2e1 * t52, t41 * t20 ^ 2 + t21 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -t40, -t8, 0, 0, t26, t17, 0, t27, 0, 0, t25 - t46, t3 - t50, t42 + t51, -t5 * pkin(4) + pkin(8) * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t24, -t47, 0, 0, t26, t17, 0, t27, 0, 0, t25 - t45, t16 - t50, t42 + t52, -t21 * pkin(4) + pkin(8) * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t26, t17, 0, t27, 0, 0, 0.2e1 * t25, -0.2e1 * t50, 0.2e1 * t42, t41 * pkin(8) ^ 2 + pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, t33, 0, -t30 * t6, -t33 * t6, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t30, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, t33, 0, -t30 * t20, -t33 * t20, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, t33, 0, -t30 * pkin(8), -t33 * pkin(8), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t1;
