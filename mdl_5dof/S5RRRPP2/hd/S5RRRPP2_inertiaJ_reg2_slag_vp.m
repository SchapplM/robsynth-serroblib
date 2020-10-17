% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPP2_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP2_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP2_inertiaJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:51:59
% EndTime: 2019-12-31 20:52:01
% DurationCPUTime: 0.62s
% Computational Cost: add. (218->61), mult. (399->91), div. (0->0), fcn. (317->4), ass. (0->51)
t38 = sin(qJ(3));
t36 = t38 ^ 2;
t40 = cos(qJ(3));
t37 = t40 ^ 2;
t67 = t36 + t37;
t39 = sin(qJ(2));
t60 = t39 * pkin(1);
t24 = pkin(7) + t60;
t53 = t67 * t24;
t66 = -0.2e1 * t38;
t65 = 0.2e1 * t38;
t64 = -0.2e1 * t40;
t63 = 0.2e1 * t40;
t33 = t40 * pkin(4);
t41 = cos(qJ(2));
t35 = t41 * pkin(1);
t48 = t40 * pkin(3) + t38 * qJ(4) + pkin(2);
t4 = -t35 - t48;
t3 = t33 - t4;
t5 = t33 + t48;
t62 = t3 + t5;
t61 = t38 * pkin(7);
t25 = -t35 - pkin(2);
t59 = pkin(2) - t25;
t58 = t48 - t4;
t57 = t38 * t24;
t56 = t38 * t40;
t55 = t53 * pkin(7);
t54 = t67 * t24 ^ 2;
t52 = t67 * pkin(7) ^ 2;
t51 = t67 * pkin(7);
t50 = t40 * qJ(4);
t49 = t40 * qJ(5);
t14 = -t38 * pkin(3) + t50;
t45 = qJ(4) ^ 2;
t44 = 0.2e1 * qJ(4);
t42 = pkin(3) + pkin(4);
t32 = t40 * pkin(7);
t26 = t38 * qJ(5);
t22 = -0.2e1 * t56;
t21 = 0.2e1 * t56;
t20 = t40 * t24;
t15 = t32 - t49;
t13 = -t26 + t61;
t9 = t42 * t38 - t50;
t8 = 0.2e1 * t51;
t7 = t20 - t49;
t6 = -t26 + t57;
t2 = 0.2e1 * t53;
t1 = t51 + t53;
t10 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t35, -0.2e1 * t60, 0, (t39 ^ 2 + t41 ^ 2) * pkin(1) ^ 2, t36, t21, 0, t37, 0, 0, t25 * t64, t25 * t65, t2, t25 ^ 2 + t54, t36, 0, t22, 0, 0, t37, t4 * t64, t2, t4 * t66, t4 ^ 2 + t54, t36, t22, 0, t37, 0, 0, t3 * t63, t3 * t65, -0.2e1 * t6 * t38 - 0.2e1 * t7 * t40, t3 ^ 2 + t6 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t35, -t60, 0, 0, t36, t21, 0, t37, 0, 0, t59 * t40, -t59 * t38, t1, -t25 * pkin(2) + t55, t36, 0, t22, 0, 0, t37, t58 * t40, t1, t58 * t38, -t4 * t48 + t55, t36, t22, 0, t37, 0, 0, t62 * t40, t62 * t38, (-t15 - t7) * t40 + (-t13 - t6) * t38, t6 * t13 + t7 * t15 + t3 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t36, t21, 0, t37, 0, 0, pkin(2) * t63, pkin(2) * t66, t8, pkin(2) ^ 2 + t52, t36, 0, t22, 0, 0, t37, -t48 * t64, t8, -t48 * t66, t48 ^ 2 + t52, t36, t22, 0, t37, 0, 0, t5 * t63, t5 * t65, -0.2e1 * t13 * t38 - 0.2e1 * t15 * t40, t13 ^ 2 + t15 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, t40, 0, -t57, -t20, 0, 0, 0, t38, 0, 0, -t40, 0, -t57, t14, t20, t14 * t24, 0, 0, -t38, 0, t40, 0, -t6, t7, t9, t7 * qJ(4) - t6 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, t40, 0, -t61, -t32, 0, 0, 0, t38, 0, 0, -t40, 0, -t61, t14, t32, t14 * pkin(7), 0, 0, -t38, 0, t40, 0, -t13, t15, t9, t15 * qJ(4) - t13 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3), 0, t44, pkin(3) ^ 2 + t45, 0, 0, 0, 0, 0, 1, 0.2e1 * t42, t44, 0, t42 ^ 2 + t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, t57, 0, 0, 0, 0, 0, 0, 0, 0, -t38, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, t61, 0, 0, 0, 0, 0, 0, 0, 0, -t38, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(3), 0, 0, 0, 0, 0, 0, -1, 0, 0, -t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t38, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t38, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t10;
