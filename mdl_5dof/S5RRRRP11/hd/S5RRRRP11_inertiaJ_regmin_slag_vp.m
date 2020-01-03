% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRP11_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP11_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP11_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t45 = cos(pkin(5));
t47 = sin(qJ(3));
t50 = cos(qJ(3));
t44 = sin(pkin(5));
t48 = sin(qJ(2));
t77 = t44 * t48;
t26 = t45 * t47 + t50 * t77;
t46 = sin(qJ(4));
t49 = cos(qJ(4));
t51 = cos(qJ(2));
t76 = t44 * t51;
t14 = t26 * t46 + t49 * t76;
t96 = -0.2e1 * t14;
t95 = -0.2e1 * t26;
t55 = -t49 * pkin(4) - t46 * qJ(5);
t31 = -pkin(3) + t55;
t94 = -0.2e1 * t31;
t93 = -0.2e1 * t47;
t92 = 0.2e1 * t50;
t60 = pkin(7) * t76;
t91 = pkin(1) * t48;
t22 = t60 + (pkin(8) + t91) * t45;
t23 = (-pkin(2) * t51 - pkin(8) * t48 - pkin(1)) * t44;
t13 = t50 * t22 + t47 * t23;
t11 = -pkin(9) * t76 + t13;
t34 = pkin(7) * t77;
t90 = pkin(1) * t51;
t21 = t34 + (-pkin(2) - t90) * t45;
t25 = -t45 * t50 + t47 * t77;
t9 = t25 * pkin(3) - t26 * pkin(9) + t21;
t4 = t49 * t11 + t46 * t9;
t89 = pkin(3) * t46;
t88 = pkin(3) * t49;
t87 = pkin(8) * t46;
t86 = pkin(8) * t49;
t85 = t25 * pkin(4);
t84 = t46 * pkin(9);
t83 = t49 * pkin(9);
t12 = -t47 * t22 + t50 * t23;
t10 = pkin(3) * t76 - t12;
t82 = t10 * t46;
t81 = t10 * t49;
t80 = t14 * t49;
t15 = t26 * t49 - t46 * t76;
t79 = t15 * t46;
t40 = t44 ^ 2;
t78 = t40 * t51;
t75 = t45 * t48;
t74 = t46 * t25;
t73 = t46 * t47;
t72 = t46 * t49;
t71 = t46 * t50;
t70 = t47 * t25;
t69 = t49 * t25;
t39 = t49 * t47;
t68 = t49 * t50;
t32 = -t50 * pkin(3) - t47 * pkin(9) - pkin(2);
t20 = pkin(8) * t68 + t46 * t32;
t41 = t46 ^ 2;
t43 = t49 ^ 2;
t67 = t41 + t43;
t66 = t25 * qJ(5);
t65 = t50 * qJ(5);
t64 = 0.2e1 * t76;
t63 = t47 * t92;
t62 = pkin(9) * t74;
t61 = pkin(9) * t69;
t59 = t47 * t76;
t58 = t50 * t76;
t57 = t46 * t11 - t49 * t9;
t1 = t66 + t4;
t2 = t57 - t85;
t56 = t1 * t49 + t2 * t46;
t54 = -pkin(4) * t46 + t49 * qJ(5);
t16 = -t65 + t20;
t30 = t49 * t32;
t17 = -t30 + (pkin(4) + t87) * t50;
t53 = t16 * t49 + t17 * t46;
t42 = t47 ^ 2;
t36 = pkin(9) * t71;
t28 = pkin(1) * t75 + t60;
t27 = t45 * t90 - t34;
t24 = (pkin(8) - t54) * t47;
t19 = -pkin(8) * t71 + t30;
t5 = t14 * pkin(4) - t15 * qJ(5) + t10;
t3 = [1, 0, 0, t40 * t48 ^ 2, 0.2e1 * t48 * t78, 0.2e1 * t44 * t75, t45 * t64, t45 ^ 2, 0.2e1 * pkin(1) * t78 + 0.2e1 * t27 * t45, -0.2e1 * t28 * t45 - 0.2e1 * t40 * t91, t26 ^ 2, t25 * t95, t76 * t95, t25 * t64, t40 * t51 ^ 2, -0.2e1 * t12 * t76 + 0.2e1 * t21 * t25, 0.2e1 * t13 * t76 + 0.2e1 * t21 * t26, t15 ^ 2, t15 * t96, 0.2e1 * t15 * t25, t25 * t96, t25 ^ 2, 0.2e1 * t10 * t14 - 0.2e1 * t25 * t57, 0.2e1 * t10 * t15 - 0.2e1 * t4 * t25, 0.2e1 * t5 * t14 - 0.2e1 * t2 * t25, -0.2e1 * t1 * t14 + 0.2e1 * t2 * t15, 0.2e1 * t1 * t25 - 0.2e1 * t5 * t15, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t77, t76, t45, t27, -t28, t26 * t47, t26 * t50 - t70, -t59, -t58, 0, -pkin(2) * t25 + pkin(8) * t59 - t21 * t50, -pkin(2) * t26 + pkin(8) * t58 + t21 * t47, t15 * t39, (-t79 - t80) * t47, -t15 * t50 + t25 * t39, t14 * t50 - t46 * t70, -t25 * t50, t19 * t25 + t57 * t50 + (pkin(8) * t14 + t82) * t47, -t20 * t25 + t4 * t50 + (pkin(8) * t15 + t81) * t47, t24 * t14 - t17 * t25 + t2 * t50 + t5 * t73, -t16 * t14 + t17 * t15 + (-t1 * t46 + t2 * t49) * t47, -t1 * t50 - t24 * t15 + t16 * t25 - t5 * t39, t1 * t16 + t2 * t17 + t5 * t24; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t42, t63, 0, 0, 0, pkin(2) * t92, pkin(2) * t93, t43 * t42, -0.2e1 * t42 * t72, t68 * t93, t46 * t63, t50 ^ 2, -0.2e1 * t19 * t50 + 0.2e1 * t42 * t87, 0.2e1 * t20 * t50 + 0.2e1 * t42 * t86, 0.2e1 * t17 * t50 + 0.2e1 * t24 * t73, 0.2e1 * (-t16 * t46 + t17 * t49) * t47, -0.2e1 * t16 * t50 - 0.2e1 * t24 * t39, t16 ^ 2 + t17 ^ 2 + t24 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t25, -t76, t12, -t13, t79, -t46 * t14 + t15 * t49, t74, t69, 0, -pkin(3) * t14 - t62 - t81, -pkin(3) * t15 - t61 + t82, t31 * t14 - t5 * t49 - t62, (t79 - t80) * pkin(9) + t56, -t31 * t15 - t5 * t46 + t61, pkin(9) * t56 + t5 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, t50, 0, -t47 * pkin(8), -t50 * pkin(8), t46 * t39, (-t41 + t43) * t47, -t71, -t68, 0, t36 + (-t86 - t89) * t47, pkin(9) * t68 + (t87 - t88) * t47, -t24 * t49 + t31 * t73 + t36, t53, -t24 * t46 + (-pkin(9) * t50 - t31 * t47) * t49, pkin(9) * t53 + t24 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t41, 0.2e1 * t72, 0, 0, 0, 0.2e1 * t88, -0.2e1 * t89, t49 * t94, 0.2e1 * t67 * pkin(9), t46 * t94, t67 * pkin(9) ^ 2 + t31 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t14, t25, -t57, -t4, -t57 + 0.2e1 * t85, -t15 * pkin(4) - t14 * qJ(5), 0.2e1 * t66 + t4, -t2 * pkin(4) + t1 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t73, -t50, t19, -t20, t30 + (-0.2e1 * pkin(4) - t87) * t50, t55 * t47, -0.2e1 * t65 + t20, -t17 * pkin(4) + t16 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, t49, 0, -t84, -t83, -t84, t54, t83, t54 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, t15, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, t39, 0, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0, t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t3;
