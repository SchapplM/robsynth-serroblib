% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPRRR6
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
% Datum: 2019-12-31 19:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRR6_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t45 = sin(qJ(4));
t71 = t45 * pkin(3);
t34 = pkin(8) + t71;
t44 = sin(qJ(5));
t38 = t44 ^ 2;
t47 = cos(qJ(5));
t40 = t47 ^ 2;
t59 = t38 + t40;
t80 = t59 * t34;
t49 = cos(qJ(3));
t42 = sin(pkin(9));
t73 = t42 * pkin(1);
t31 = pkin(6) + t73;
t66 = pkin(7) + t31;
t18 = t66 * t49;
t48 = cos(qJ(4));
t46 = sin(qJ(3));
t56 = t66 * t46;
t6 = t45 * t18 + t48 * t56;
t79 = t6 ^ 2;
t22 = t45 * t46 - t48 * t49;
t78 = t22 ^ 2;
t43 = cos(pkin(9));
t72 = t43 * pkin(1);
t32 = -pkin(2) - t72;
t25 = -t49 * pkin(3) + t32;
t77 = 0.2e1 * t25;
t76 = 0.2e1 * t46;
t75 = t22 * pkin(4);
t24 = t45 * t49 + t48 * t46;
t74 = t24 * pkin(8);
t70 = t48 * pkin(3);
t69 = t6 * t22;
t68 = t6 * t47;
t35 = -pkin(4) - t70;
t67 = pkin(4) - t35;
t65 = t38 * t24;
t64 = t44 * t24;
t63 = t44 * t47;
t62 = t47 * t24;
t60 = pkin(8) * t59;
t39 = t46 ^ 2;
t41 = t49 ^ 2;
t58 = t39 + t41;
t57 = -0.2e1 * t24 * t22;
t54 = -pkin(4) * t24 - pkin(8) * t22;
t5 = t25 - t74 + t75;
t8 = t48 * t18 - t45 * t56;
t2 = -t44 * t8 + t47 * t5;
t3 = t44 * t5 + t47 * t8;
t1 = -t2 * t44 + t3 * t47;
t53 = -t22 * t34 + t24 * t35;
t29 = 0.2e1 * t63;
t21 = t24 ^ 2;
t17 = t47 * t22;
t16 = t40 * t24;
t15 = t40 * t21;
t14 = t44 * t22;
t13 = t38 * t21;
t11 = t44 * t62;
t10 = t16 + t65;
t9 = t16 - t65;
t4 = t6 * t44;
t7 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t72, -0.2e1 * t73, 0, (t42 ^ 2 + t43 ^ 2) * pkin(1) ^ 2, t39, t49 * t76, 0, t41, 0, 0, -0.2e1 * t32 * t49, t32 * t76, 0.2e1 * t58 * t31, t58 * t31 ^ 2 + t32 ^ 2, t21, t57, 0, t78, 0, 0, t22 * t77, t24 * t77, -0.2e1 * t8 * t22 + 0.2e1 * t6 * t24, t25 ^ 2 + t8 ^ 2 + t79, t15, -0.2e1 * t21 * t63, 0.2e1 * t22 * t62, t13, t44 * t57, t78, 0.2e1 * t2 * t22 + 0.2e1 * t6 * t64, -0.2e1 * t3 * t22 + 0.2e1 * t6 * t62, 0.2e1 * (-t2 * t47 - t3 * t44) * t24, t2 ^ 2 + t3 ^ 2 + t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8 * t24 + t69, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 * t24 + t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21 + t78, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15 + t13 + t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0, t49, 0, -t46 * t31, -t49 * t31, 0, 0, 0, 0, t24, 0, -t22, 0, -t6, -t8, (-t22 * t45 - t24 * t48) * pkin(3), (t45 * t8 - t48 * t6) * pkin(3), t11, t9, t14, -t11, t17, 0, t53 * t44 - t68, t53 * t47 + t4, t1, t1 * t34 + t6 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t46, 0, 0, 0, 0, 0, 0, 0, 0, -t22, -t24, 0, (-t22 * t48 + t24 * t45) * pkin(3), 0, 0, 0, 0, 0, 0, -t17, t14, t10, t22 * t35 + t24 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t70, -0.2e1 * t71, 0, (t45 ^ 2 + t48 ^ 2) * pkin(3) ^ 2, t38, t29, 0, t40, 0, 0, -0.2e1 * t35 * t47, 0.2e1 * t35 * t44, 0.2e1 * t80, t59 * t34 ^ 2 + t35 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, -t22, 0, -t6, -t8, 0, 0, t11, t9, t14, -t11, t17, 0, t54 * t44 - t68, t54 * t47 + t4, t1, -t6 * pkin(4) + t1 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, -t24, 0, 0, 0, 0, 0, 0, 0, 0, -t17, t14, t10, t59 * t74 - t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t70, -t71, 0, 0, t38, t29, 0, t40, 0, 0, t67 * t47, -t67 * t44, t60 + t80, -t35 * pkin(4) + pkin(8) * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t38, t29, 0, t40, 0, 0, 0.2e1 * pkin(4) * t47, -0.2e1 * pkin(4) * t44, 0.2e1 * t60, t59 * pkin(8) ^ 2 + pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, -t64, t22, t2, -t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, -t62, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, t47, 0, -t44 * t34, -t47 * t34, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, t47, 0, -t44 * pkin(8), -t47 * pkin(8), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t7;
