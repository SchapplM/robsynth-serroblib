% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:52
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPPR1_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:51:43
% EndTime: 2022-01-20 09:51:46
% DurationCPUTime: 0.60s
% Computational Cost: add. (379->56), mult. (680->99), div. (0->0), fcn. (682->8), ass. (0->52)
t44 = sin(pkin(8));
t63 = t44 * pkin(2);
t35 = qJ(4) + t63;
t43 = sin(pkin(9));
t41 = t43 ^ 2;
t45 = cos(pkin(9));
t42 = t45 ^ 2;
t54 = t41 + t42;
t56 = t54 * t35;
t52 = -t45 * pkin(4) - pkin(3);
t49 = cos(qJ(2));
t40 = t49 * pkin(1);
t38 = t40 + pkin(2);
t46 = cos(pkin(8));
t48 = sin(qJ(2));
t61 = t48 * pkin(1);
t55 = -t46 * t38 + t44 * t61;
t13 = t52 + t55;
t69 = 0.2e1 * t13;
t62 = t46 * pkin(2);
t28 = t52 - t62;
t68 = 0.2e1 * t28;
t67 = 0.2e1 * t43;
t66 = -0.2e1 * t45;
t47 = sin(qJ(5));
t60 = cos(qJ(5));
t25 = t47 * t43 - t60 * t45;
t27 = t60 * t43 + t47 * t45;
t53 = t46 * t61;
t20 = t44 * t38 + t53;
t17 = qJ(4) + t20;
t11 = (-pkin(7) - t17) * t43;
t39 = t45 * pkin(7);
t12 = t45 * t17 + t39;
t3 = t60 * t11 - t47 * t12;
t4 = t47 * t11 + t60 * t12;
t65 = -t4 * t25 - t3 * t27;
t21 = (-pkin(7) - t35) * t43;
t22 = t45 * t35 + t39;
t8 = t60 * t21 - t47 * t22;
t9 = t47 * t21 + t60 * t22;
t64 = -t9 * t25 - t8 * t27;
t59 = t13 + t28;
t58 = t54 * t17;
t18 = -pkin(3) + t55;
t37 = -pkin(3) - t62;
t57 = t18 + t37;
t32 = t45 * t67;
t24 = t27 ^ 2;
t23 = t25 ^ 2;
t10 = -0.2e1 * t27 * t25;
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t40, -0.2e1 * t61, 0, (t48 ^ 2 + t49 ^ 2) * pkin(1) ^ 2, 0, 0, 0, 0, 0, 1, -0.2e1 * t55, -0.2e1 * t20, 0, t20 ^ 2 + t55 ^ 2, t41, t32, 0, t42, 0, 0, t18 * t66, t18 * t67, 0.2e1 * t58, t54 * t17 ^ 2 + t18 ^ 2, t24, t10, 0, t23, 0, 0, t25 * t69, t27 * t69, 0.2e1 * t65, t13 ^ 2 + t3 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t40, -t61, 0, 0, 0, 0, 0, 0, 0, 1, -t55 + t62, -t53 + (-pkin(2) - t38) * t44, 0, (t20 * t44 - t46 * t55) * pkin(2), t41, t32, 0, t42, 0, 0, -t57 * t45, t57 * t43, t56 + t58, t17 * t56 + t18 * t37, t24, t10, 0, t23, 0, 0, t59 * t25, t59 * t27, t64 + t65, t13 * t28 + t3 * t8 + t4 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t62, -0.2e1 * t63, 0, (t44 ^ 2 + t46 ^ 2) * pkin(2) ^ 2, t41, t32, 0, t42, 0, 0, t37 * t66, t37 * t67, 0.2e1 * t56, t54 * t35 ^ 2 + t37 ^ 2, t24, t10, 0, t23, 0, 0, t25 * t68, t27 * t68, 0.2e1 * t64, t28 ^ 2 + t8 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3 * t25 + t4 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8 * t25 + t9 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24 + t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, t43, 0, t18, 0, 0, 0, 0, 0, 0, t25, t27, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, t43, 0, t37, 0, 0, 0, 0, 0, 0, t25, t27, 0, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, -t25, 0, t3, -t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, -t25, 0, t8, -t9, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, -t27, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t1;
