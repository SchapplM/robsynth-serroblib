% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRP9_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP9_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP9_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t47 = sin(pkin(8));
t48 = cos(pkin(8));
t49 = sin(qJ(4));
t82 = cos(qJ(4));
t32 = t82 * t47 + t49 * t48;
t40 = -t48 * pkin(3) - pkin(2);
t88 = -t49 * t47 + t82 * t48;
t10 = -pkin(4) * t88 - t32 * qJ(5) + t40;
t89 = -0.2e1 * t10;
t50 = sin(qJ(2));
t23 = t32 * t50;
t87 = t23 ^ 2;
t86 = t88 ^ 2;
t85 = 0.2e1 * t40;
t45 = t50 ^ 2;
t84 = t45 * pkin(6);
t42 = t50 * pkin(6);
t51 = cos(qJ(2));
t83 = t51 * pkin(4);
t34 = -t51 * pkin(2) - t50 * qJ(3) - pkin(1);
t27 = t48 * t34;
t72 = t48 * t50;
t11 = -pkin(7) * t72 + t27 + (-pkin(6) * t47 - pkin(3)) * t51;
t71 = t48 * t51;
t21 = pkin(6) * t71 + t47 * t34;
t74 = t47 * t50;
t14 = -pkin(7) * t74 + t21;
t4 = t49 * t11 + t82 * t14;
t66 = pkin(7) + qJ(3);
t35 = t66 * t48;
t58 = t66 * t47;
t16 = t49 * t35 + t82 * t58;
t81 = t16 * t51;
t18 = t82 * t35 - t49 * t58;
t80 = t18 * t51;
t79 = t23 * t88;
t25 = t88 * t50;
t78 = t25 * t23;
t77 = t32 * t88;
t76 = t32 * t51;
t75 = t47 * t48;
t73 = t47 * t51;
t69 = t50 * t51;
t68 = t51 * t23;
t67 = t51 * t88;
t33 = pkin(3) * t74 + t42;
t43 = t47 ^ 2;
t44 = t48 ^ 2;
t65 = t43 + t44;
t64 = t51 * qJ(5);
t63 = 0.2e1 * t69;
t62 = t16 ^ 2 + t18 ^ 2;
t61 = t47 * t72;
t3 = t82 * t11 - t49 * t14;
t59 = t16 * t25 - t18 * t23;
t57 = -pkin(2) * t50 + qJ(3) * t51;
t20 = -pkin(6) * t73 + t27;
t56 = -t20 * t47 + t21 * t48;
t55 = -t32 * t23 + t25 * t88;
t54 = 0.2e1 * t16 * t32 + 0.2e1 * t18 * t88;
t53 = pkin(6) ^ 2;
t46 = t51 ^ 2;
t41 = t45 * t53;
t28 = t32 ^ 2;
t22 = t25 ^ 2;
t19 = -0.2e1 * t25 * t51;
t12 = t25 * t32;
t5 = t23 * pkin(4) - t25 * qJ(5) + t33;
t2 = -t3 + t83;
t1 = -t64 + t4;
t6 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t45, t63, 0, t46, 0, 0, 0.2e1 * pkin(1) * t51, -0.2e1 * pkin(1) * t50, 0.2e1 * (t45 + t46) * pkin(6), pkin(1) ^ 2 + t46 * t53 + t41, t44 * t45, -0.2e1 * t45 * t75, -0.2e1 * t48 * t69, t43 * t45, t47 * t63, t46, -0.2e1 * t20 * t51 + 0.2e1 * t47 * t84, 0.2e1 * t21 * t51 + 0.2e1 * t48 * t84, 0.2e1 * (-t20 * t48 - t21 * t47) * t50, t20 ^ 2 + t21 ^ 2 + t41, t22, -0.2e1 * t78, t19, t87, 0.2e1 * t68, t46, 0.2e1 * t33 * t23 - 0.2e1 * t3 * t51, 0.2e1 * t33 * t25 + 0.2e1 * t4 * t51, -0.2e1 * t4 * t23 - 0.2e1 * t3 * t25, t3 ^ 2 + t33 ^ 2 + t4 ^ 2, t22, t19, 0.2e1 * t78, t46, -0.2e1 * t68, t87, 0.2e1 * t2 * t51 + 0.2e1 * t5 * t23, -0.2e1 * t1 * t23 + 0.2e1 * t2 * t25, -0.2e1 * t1 * t51 - 0.2e1 * t5 * t25, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, t51, 0, -t42, -t51 * pkin(6), 0, 0, t61, (-t43 + t44) * t50, -t73, -t61, -t71, 0, -pkin(6) * t72 + t57 * t47, pkin(6) * t74 + t57 * t48, t56, -pkin(2) * t42 + t56 * qJ(3), t12, t55, -t76, -t79, -t67, 0, t40 * t23 - t33 * t88 + t81, t40 * t25 + t33 * t32 + t80, -t3 * t32 + t4 * t88 + t59, -t3 * t16 + t4 * t18 + t33 * t40, t12, -t76, -t55, 0, t67, -t79, t10 * t23 - t5 * t88 + t81, t1 * t88 + t2 * t32 + t59, -t10 * t25 - t5 * t32 - t80, t1 * t18 + t5 * t10 + t2 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t43, 0.2e1 * t75, 0, t44, 0, 0, 0.2e1 * pkin(2) * t48, -0.2e1 * pkin(2) * t47, 0.2e1 * t65 * qJ(3), t65 * qJ(3) ^ 2 + pkin(2) ^ 2, t28, 0.2e1 * t77, 0, t86, 0, 0, -t88 * t85, t32 * t85, t54, t40 ^ 2 + t62, t28, 0, -0.2e1 * t77, 0, 0, t86, t88 * t89, t54, t32 * t89, t10 ^ 2 + t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, t72, 0, t42, 0, 0, 0, 0, 0, 0, t23, t25, 0, t33, 0, 0, 0, 0, 0, 0, t23, 0, -t25, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, t47, 0, -pkin(2), 0, 0, 0, 0, 0, 0, -t88, t32, 0, t40, 0, 0, 0, 0, 0, 0, -t88, 0, -t32, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, -t23, -t51, t3, -t4, 0, 0, 0, t25, 0, -t51, t23, 0, t3 - 0.2e1 * t83, -pkin(4) * t25 - t23 * qJ(5), -0.2e1 * t64 + t4, -t2 * pkin(4) + t1 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, t88, 0, -t16, -t18, 0, 0, 0, t32, 0, 0, -t88, 0, -t16, -pkin(4) * t32 + qJ(5) * t88, t18, -t16 * pkin(4) + t18 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, t25, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t6;
