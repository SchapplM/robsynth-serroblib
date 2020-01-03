% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRP9_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP9_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP9_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t49 = sin(qJ(4));
t50 = sin(qJ(3));
t52 = cos(qJ(4));
t53 = cos(qJ(3));
t30 = t49 * t53 + t52 * t50;
t51 = sin(qJ(2));
t23 = t30 * t51;
t96 = t23 ^ 2;
t28 = t49 * t50 - t52 * t53;
t95 = t28 ^ 2;
t41 = -t53 * pkin(3) - pkin(2);
t94 = 0.2e1 * t41;
t93 = -0.2e1 * t51;
t54 = cos(qJ(2));
t92 = 0.2e1 * t54;
t91 = -pkin(8) - pkin(7);
t90 = pkin(2) * t53;
t89 = pkin(6) * t50;
t46 = t51 ^ 2;
t88 = t46 * pkin(6);
t44 = t51 * pkin(6);
t87 = t52 * pkin(3);
t86 = t54 * pkin(3);
t85 = t54 * pkin(4);
t33 = -t54 * pkin(2) - t51 * pkin(7) - pkin(1);
t26 = t53 * t33;
t75 = t53 * t51;
t11 = -pkin(8) * t75 + t26 + (-pkin(3) - t89) * t54;
t74 = t53 * t54;
t67 = pkin(6) * t74;
t14 = t67 + (-pkin(8) * t51 + t33) * t50;
t4 = t49 * t11 + t52 * t14;
t34 = t91 * t53;
t65 = t91 * t50;
t16 = -t49 * t34 - t52 * t65;
t84 = t16 * t54;
t18 = -t52 * t34 + t49 * t65;
t83 = t18 * t54;
t82 = t23 * t28;
t78 = t50 * t51;
t25 = -t49 * t78 + t52 * t75;
t81 = t25 * t23;
t80 = t30 * t28;
t79 = t30 * t54;
t77 = t50 * t53;
t76 = t50 * t54;
t73 = t54 * t23;
t72 = t54 * t28;
t32 = pkin(3) * t78 + t44;
t45 = t50 ^ 2;
t47 = t53 ^ 2;
t71 = t45 + t47;
t70 = t54 * qJ(5);
t69 = t51 * t92;
t68 = t16 ^ 2 + t18 ^ 2;
t66 = t50 * t75;
t64 = -t52 * t11 + t49 * t14;
t63 = t16 * t25 - t18 * t23;
t20 = -pkin(6) * t76 + t26;
t21 = t50 * t33 + t67;
t62 = -t20 * t50 + t21 * t53;
t61 = -t30 * t23 - t25 * t28;
t60 = 0.2e1 * t16 * t30 - 0.2e1 * t18 * t28;
t58 = pkin(6) ^ 2;
t56 = 0.2e1 * pkin(4);
t55 = 0.2e1 * qJ(5);
t48 = t54 ^ 2;
t43 = t46 * t58;
t42 = t49 * pkin(3);
t39 = pkin(4) + t87;
t37 = t42 + qJ(5);
t27 = t30 ^ 2;
t22 = t25 ^ 2;
t19 = -0.2e1 * t25 * t54;
t12 = t25 * t30;
t10 = t28 * pkin(4) - t30 * qJ(5) + t41;
t5 = t23 * pkin(4) - t25 * qJ(5) + t32;
t2 = t64 + t85;
t1 = -t70 + t4;
t3 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t46, t69, 0, t48, 0, 0, pkin(1) * t92, pkin(1) * t93, 0.2e1 * (t46 + t48) * pkin(6), pkin(1) ^ 2 + t48 * t58 + t43, t47 * t46, -0.2e1 * t46 * t77, t74 * t93, t45 * t46, t50 * t69, t48, -0.2e1 * t20 * t54 + 0.2e1 * t50 * t88, 0.2e1 * t21 * t54 + 0.2e1 * t53 * t88, 0.2e1 * (-t20 * t53 - t21 * t50) * t51, t20 ^ 2 + t21 ^ 2 + t43, t22, -0.2e1 * t81, t19, t96, 0.2e1 * t73, t48, 0.2e1 * t32 * t23 + 0.2e1 * t54 * t64, 0.2e1 * t32 * t25 + 0.2e1 * t4 * t54, -0.2e1 * t4 * t23 + 0.2e1 * t25 * t64, t32 ^ 2 + t4 ^ 2 + t64 ^ 2, t22, t19, 0.2e1 * t81, t48, -0.2e1 * t73, t96, 0.2e1 * t2 * t54 + 0.2e1 * t5 * t23, -0.2e1 * t1 * t23 + 0.2e1 * t2 * t25, -0.2e1 * t1 * t54 - 0.2e1 * t5 * t25, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, t54, 0, -t44, -t54 * pkin(6), 0, 0, t66, (-t45 + t47) * t51, -t76, -t66, -t74, 0, -pkin(6) * t75 + (-pkin(2) * t51 + pkin(7) * t54) * t50, pkin(7) * t74 + (t89 - t90) * t51, t62, -pkin(2) * t44 + pkin(7) * t62, t12, t61, -t79, t82, t72, 0, t41 * t23 + t32 * t28 + t84, t41 * t25 + t32 * t30 + t83, -t4 * t28 + t30 * t64 + t63, t16 * t64 + t4 * t18 + t32 * t41, t12, -t79, -t61, 0, -t72, t82, t10 * t23 + t5 * t28 + t84, -t1 * t28 + t2 * t30 + t63, -t10 * t25 - t5 * t30 - t83, t1 * t18 + t5 * t10 + t2 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t45, 0.2e1 * t77, 0, t47, 0, 0, 0.2e1 * t90, -0.2e1 * pkin(2) * t50, 0.2e1 * t71 * pkin(7), t71 * pkin(7) ^ 2 + pkin(2) ^ 2, t27, -0.2e1 * t80, 0, t95, 0, 0, t28 * t94, t30 * t94, t60, t41 ^ 2 + t68, t27, 0, 0.2e1 * t80, 0, 0, t95, 0.2e1 * t10 * t28, t60, -0.2e1 * t10 * t30, t10 ^ 2 + t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, 0, -t78, -t54, t20, -t21, 0, 0, 0, 0, t25, 0, -t23, -t54, -t52 * t86 - t64, t49 * t86 - t4, (-t23 * t49 - t25 * t52) * pkin(3), (t4 * t49 - t52 * t64) * pkin(3), 0, t25, 0, -t54, t23, 0, (-pkin(4) - t39) * t54 - t64, -t37 * t23 - t39 * t25, (-qJ(5) - t37) * t54 + t4, t1 * t37 - t2 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, t53, 0, -t50 * pkin(7), -t53 * pkin(7), 0, 0, 0, 0, t30, 0, -t28, 0, -t16, -t18, (-t28 * t49 - t30 * t52) * pkin(3), (-t16 * t52 + t18 * t49) * pkin(3), 0, t30, 0, 0, t28, 0, -t16, -t37 * t28 - t39 * t30, t18, -t16 * t39 + t18 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t87, -0.2e1 * t42, 0, (t49 ^ 2 + t52 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 1, 0, 0, 0.2e1 * t39, 0, 0.2e1 * t37, t37 ^ 2 + t39 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, -t23, -t54, -t64, -t4, 0, 0, 0, t25, 0, -t54, t23, 0, -t64 - 0.2e1 * t85, -pkin(4) * t25 - t23 * qJ(5), -0.2e1 * t70 + t4, -t2 * pkin(4) + t1 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, -t28, 0, -t16, -t18, 0, 0, 0, t30, 0, 0, t28, 0, -t16, -pkin(4) * t30 - t28 * qJ(5), t18, -t16 * pkin(4) + t18 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t87, -t42, 0, 0, 0, 0, 0, 1, 0, 0, t56 + t87, 0, t55 + t42, t39 * pkin(4) + t37 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t56, 0, t55, pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, t25, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t3;
