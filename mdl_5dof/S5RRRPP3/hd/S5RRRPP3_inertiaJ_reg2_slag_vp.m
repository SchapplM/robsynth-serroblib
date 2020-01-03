% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRRPP3
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
% Datum: 2019-12-31 20:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPP3_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_inertiaJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t43 = sin(qJ(3));
t39 = t43 ^ 2;
t45 = cos(qJ(3));
t40 = t45 ^ 2;
t72 = t39 + t40;
t44 = sin(qJ(2));
t63 = t44 * pkin(1);
t29 = pkin(7) + t63;
t55 = t72 * t29;
t31 = t45 * qJ(4);
t71 = -t43 * pkin(3) + t31;
t70 = -0.2e1 * t43;
t69 = -0.2e1 * t45;
t68 = 0.2e1 * t45;
t41 = pkin(3) + qJ(5);
t52 = -t43 * qJ(4) - pkin(2);
t6 = -t41 * t45 + t52;
t46 = cos(qJ(2));
t62 = t46 * pkin(1);
t3 = t6 - t62;
t67 = -t3 - t6;
t22 = t43 * t29;
t64 = t43 * pkin(4);
t8 = t22 + t64;
t25 = t45 * t29;
t38 = t45 * pkin(4);
t9 = t25 + t38;
t66 = t8 * t43 + t9 * t45;
t30 = -pkin(2) - t62;
t61 = pkin(2) - t30;
t14 = -t45 * pkin(3) + t52;
t7 = t14 - t62;
t60 = t14 + t7;
t59 = t43 * t45;
t34 = t43 * pkin(7);
t18 = t34 + t64;
t37 = t45 * pkin(7);
t19 = t37 + t38;
t58 = t18 * t43 + t19 * t45;
t57 = t55 * pkin(7);
t56 = t72 * t29 ^ 2;
t54 = t72 * pkin(7) ^ 2;
t53 = t72 * pkin(7);
t48 = qJ(4) ^ 2;
t47 = 0.2e1 * qJ(4);
t27 = -0.2e1 * t59;
t26 = 0.2e1 * t59;
t13 = -t41 * t43 + t31;
t12 = 0.2e1 * t53;
t2 = 0.2e1 * t55;
t1 = t53 + t55;
t4 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t62, -0.2e1 * t63, 0, (t44 ^ 2 + t46 ^ 2) * pkin(1) ^ 2, t39, t26, 0, t40, 0, 0, t30 * t69, 0.2e1 * t30 * t43, t2, t30 ^ 2 + t56, 0, 0, 0, t39, t26, t40, t2, t7 * t68, t7 * t70, t7 ^ 2 + t56, 0, 0, 0, t40, t27, t39, 0.2e1 * t66, t3 * t70, t3 * t69, t3 ^ 2 + t8 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t62, -t63, 0, 0, t39, t26, 0, t40, 0, 0, t61 * t45, -t61 * t43, t1, -t30 * pkin(2) + t57, 0, 0, 0, t39, t26, t40, t1, t60 * t45, -t60 * t43, t7 * t14 + t57, 0, 0, 0, t40, t27, t39, t58 + t66, t67 * t43, t67 * t45, t8 * t18 + t9 * t19 + t3 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t39, t26, 0, t40, 0, 0, pkin(2) * t68, pkin(2) * t70, t12, pkin(2) ^ 2 + t54, 0, 0, 0, t39, t26, t40, t12, t14 * t68, t14 * t70, t14 ^ 2 + t54, 0, 0, 0, t40, t27, t39, 0.2e1 * t58, t6 * t70, t6 * t69, t18 ^ 2 + t19 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, 0, t45, 0, -t22, -t25, 0, 0, 0, -t43, -t45, 0, 0, 0, t71, t22, t25, t71 * t29, 0, -t45, t43, 0, 0, 0, t13, t9, -t8, t9 * qJ(4) - t8 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, 0, t45, 0, -t34, -t37, 0, 0, 0, -t43, -t45, 0, 0, 0, t71, t34, t37, t71 * pkin(7), 0, -t45, t43, 0, 0, 0, t13, t19, -t18, t19 * qJ(4) - t18 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -0.2e1 * pkin(3), t47, pkin(3) ^ 2 + t48, 1, 0, 0, 0, 0, 0, 0, t47, 0.2e1 * t41, t41 ^ 2 + t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, 0, 0, t22, 0, 0, 0, 0, 0, 0, t43, 0, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, 0, 0, t34, 0, 0, 0, 0, 0, 0, t43, 0, 0, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(3), 0, 0, 0, 0, 0, 0, 0, 0, -1, -t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, 0, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, qJ(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t4;
