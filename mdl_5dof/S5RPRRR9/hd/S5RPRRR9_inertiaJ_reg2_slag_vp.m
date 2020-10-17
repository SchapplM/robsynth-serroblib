% Calculate inertial parameters regressor of joint inertia matrix for
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
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRR9_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR9_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR9_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:08:00
% EndTime: 2019-12-31 19:08:02
% DurationCPUTime: 0.74s
% Computational Cost: add. (923->74), mult. (1787->155), div. (0->0), fcn. (2123->8), ass. (0->62)
t42 = sin(qJ(4));
t68 = t42 * pkin(3);
t31 = pkin(8) + t68;
t41 = sin(qJ(5));
t37 = t41 ^ 2;
t44 = cos(qJ(5));
t38 = t44 ^ 2;
t57 = t37 + t38;
t60 = t57 * t31;
t39 = sin(pkin(9));
t40 = cos(pkin(9));
t43 = sin(qJ(3));
t46 = cos(qJ(3));
t23 = -t43 * t39 + t46 * t40;
t24 = t46 * t39 + t43 * t40;
t45 = cos(qJ(4));
t18 = t42 * t23 + t45 * t24;
t74 = -0.2e1 * t18;
t61 = pkin(6) + qJ(2);
t54 = t61 * t40;
t55 = t61 * t39;
t20 = -t43 * t55 + t46 * t54;
t11 = t23 * pkin(7) + t20;
t19 = -t43 * t54 - t46 * t55;
t50 = -t24 * pkin(7) + t19;
t5 = t42 * t11 - t45 * t50;
t73 = t5 ^ 2;
t16 = -t45 * t23 + t42 * t24;
t72 = t16 ^ 2;
t30 = -t40 * pkin(2) - pkin(1);
t21 = -t23 * pkin(3) + t30;
t71 = 0.2e1 * t21;
t70 = 0.2e1 * t24;
t69 = 0.2e1 * t40;
t67 = t45 * pkin(3);
t66 = t5 * t44;
t32 = -pkin(4) - t67;
t65 = pkin(4) - t32;
t13 = t41 * t16;
t64 = t41 * t18;
t63 = t41 * t44;
t62 = t44 * t18;
t59 = t57 * pkin(8);
t35 = t39 ^ 2;
t36 = t40 ^ 2;
t58 = t35 + t36;
t56 = t16 * t74;
t53 = -pkin(4) * t18 - pkin(8) * t16;
t7 = t45 * t11 + t42 * t50;
t8 = t16 * pkin(4) - t18 * pkin(8) + t21;
t2 = -t41 * t7 + t44 * t8;
t3 = t41 * t8 + t44 * t7;
t52 = t2 * t44 + t3 * t41;
t1 = -t2 * t41 + t3 * t44;
t51 = -t16 * t31 + t18 * t32;
t27 = 0.2e1 * t63;
t15 = t18 ^ 2;
t14 = t44 * t16;
t12 = t41 * t62;
t9 = (-t37 + t38) * t18;
t4 = t5 * t41;
t6 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t35, t39 * t69, 0, t36, 0, 0, pkin(1) * t69, -0.2e1 * pkin(1) * t39, 0.2e1 * t58 * qJ(2), t58 * qJ(2) ^ 2 + pkin(1) ^ 2, t24 ^ 2, t23 * t70, 0, t23 ^ 2, 0, 0, -0.2e1 * t30 * t23, t30 * t70, -0.2e1 * t19 * t24 + 0.2e1 * t20 * t23, t19 ^ 2 + t20 ^ 2 + t30 ^ 2, t15, t56, 0, t72, 0, 0, t16 * t71, t18 * t71, -0.2e1 * t7 * t16 + 0.2e1 * t5 * t18, t21 ^ 2 + t7 ^ 2 + t73, t38 * t15, -0.2e1 * t15 * t63, 0.2e1 * t16 * t62, t37 * t15, t41 * t56, t72, 0.2e1 * t2 * t16 + 0.2e1 * t5 * t64, -0.2e1 * t3 * t16 + 0.2e1 * t5 * t62, t52 * t74, t2 ^ 2 + t3 ^ 2 + t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, t39, 0, -pkin(1), 0, 0, 0, 0, 0, 0, -t23, t24, 0, t30, 0, 0, 0, 0, 0, 0, t16, t18, 0, t21, 0, 0, 0, 0, 0, 0, t14, -t13, -t57 * t18, t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, t23, 0, t19, -t20, 0, 0, 0, 0, t18, 0, -t16, 0, -t5, -t7, (-t16 * t42 - t18 * t45) * pkin(3), (t42 * t7 - t45 * t5) * pkin(3), t12, t9, t13, -t12, t14, 0, t51 * t41 - t66, t51 * t44 + t4, t1, t1 * t31 + t5 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t67, -0.2e1 * t68, 0, (t42 ^ 2 + t45 ^ 2) * pkin(3) ^ 2, t37, t27, 0, t38, 0, 0, -0.2e1 * t32 * t44, 0.2e1 * t32 * t41, 0.2e1 * t60, t57 * t31 ^ 2 + t32 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, -t16, 0, -t5, -t7, 0, 0, t12, t9, t13, -t12, t14, 0, t53 * t41 - t66, t53 * t44 + t4, t1, -t5 * pkin(4) + t1 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t67, -t68, 0, 0, t37, t27, 0, t38, 0, 0, t65 * t44, -t65 * t41, t59 + t60, -t32 * pkin(4) + pkin(8) * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t37, t27, 0, t38, 0, 0, 0.2e1 * pkin(4) * t44, -0.2e1 * pkin(4) * t41, 0.2e1 * t59, t57 * pkin(8) ^ 2 + pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, -t64, t16, t2, -t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t41, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, t44, 0, -t41 * t31, -t44 * t31, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, t44, 0, -t41 * pkin(8), -t44 * pkin(8), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t6;
