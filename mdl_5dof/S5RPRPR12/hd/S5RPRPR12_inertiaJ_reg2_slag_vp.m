% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRPR12_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR12_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR12_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t47 = sin(pkin(8));
t49 = cos(pkin(8));
t51 = sin(qJ(3));
t73 = cos(qJ(3));
t34 = t73 * t47 + t51 * t49;
t80 = -0.2e1 * t34;
t46 = sin(pkin(9));
t48 = cos(pkin(9));
t50 = sin(qJ(5));
t72 = cos(qJ(5));
t79 = -t50 * t46 + t72 * t48;
t64 = pkin(6) + qJ(2);
t37 = t64 * t49;
t57 = t64 * t47;
t17 = t51 * t37 + t73 * t57;
t78 = t17 ^ 2;
t31 = t51 * t47 - t73 * t49;
t25 = t31 ^ 2;
t77 = 0.2e1 * t31;
t40 = -t48 * pkin(4) - pkin(3);
t76 = 0.2e1 * t40;
t41 = -t49 * pkin(2) - pkin(1);
t75 = 0.2e1 * t41;
t74 = 0.2e1 * t49;
t11 = t79 * t34;
t71 = t11 * t79;
t33 = t72 * t46 + t50 * t48;
t70 = t33 * t31;
t69 = t46 * t31;
t68 = t46 * t34;
t67 = t46 * t48;
t66 = t48 * t34;
t63 = pkin(7) + qJ(4);
t14 = t31 * pkin(3) - t34 * qJ(4) + t41;
t20 = t73 * t37 - t51 * t57;
t6 = t46 * t14 + t48 * t20;
t42 = t46 ^ 2;
t44 = t48 ^ 2;
t62 = t42 + t44;
t43 = t47 ^ 2;
t45 = t49 ^ 2;
t61 = t43 + t45;
t60 = t31 * t80;
t59 = t46 * t66;
t5 = t48 * t14 - t46 * t20;
t56 = t6 * t46 + t5 * t48;
t55 = -t5 * t46 + t6 * t48;
t54 = -pkin(3) * t34 - qJ(4) * t31;
t36 = t63 * t48;
t35 = t63 * t46;
t27 = t34 ^ 2;
t26 = t33 ^ 2;
t24 = t79 ^ 2;
t23 = t48 * t31;
t21 = t79 * t31;
t19 = -t50 * t35 + t72 * t36;
t16 = -t72 * t35 - t50 * t36;
t9 = t33 * t34;
t8 = pkin(4) * t68 + t17;
t7 = t33 * t9;
t4 = -pkin(7) * t68 + t6;
t3 = t31 * pkin(4) - pkin(7) * t66 + t5;
t2 = t50 * t3 + t72 * t4;
t1 = t72 * t3 - t50 * t4;
t10 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t43, t47 * t74, 0, t45, 0, 0, pkin(1) * t74, -0.2e1 * pkin(1) * t47, 0.2e1 * t61 * qJ(2), t61 * qJ(2) ^ 2 + pkin(1) ^ 2, t27, t60, 0, t25, 0, 0, t31 * t75, t34 * t75, 0.2e1 * t17 * t34 - 0.2e1 * t20 * t31, t20 ^ 2 + t41 ^ 2 + t78, t44 * t27, -0.2e1 * t27 * t67, t66 * t77, t42 * t27, t46 * t60, t25, 0.2e1 * t17 * t68 + 0.2e1 * t5 * t31, 0.2e1 * t17 * t66 - 0.2e1 * t6 * t31, t56 * t80, t5 ^ 2 + t6 ^ 2 + t78, t11 ^ 2, -0.2e1 * t11 * t9, t11 * t77, t9 ^ 2, -t9 * t77, t25, 0.2e1 * t1 * t31 + 0.2e1 * t8 * t9, 0.2e1 * t8 * t11 - 0.2e1 * t2 * t31, -0.2e1 * t1 * t11 - 0.2e1 * t2 * t9, t1 ^ 2 + t2 ^ 2 + t8 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, t47, 0, -pkin(1), 0, 0, 0, 0, 0, 0, t31, t34, 0, t41, 0, 0, 0, 0, 0, 0, t23, -t69, -t62 * t34, t56, 0, 0, 0, 0, 0, 0, t21, -t70, -t7 - t71, t1 * t79 + t2 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26 + t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, -t31, 0, -t17, -t20, 0, 0, t59, (-t42 + t44) * t34, t69, -t59, t23, 0, -t17 * t48 + t54 * t46, t17 * t46 + t54 * t48, t55, -t17 * pkin(3) + t55 * qJ(4), t11 * t33, -t7 + t71, t70, -t9 * t79, t21, 0, t16 * t31 + t40 * t9 - t79 * t8, t40 * t11 - t19 * t31 + t8 * t33, -t1 * t33 - t16 * t11 - t19 * t9 + t2 * t79, t1 * t16 + t2 * t19 + t8 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16 * t79 + t33 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t42, 0.2e1 * t67, 0, t44, 0, 0, 0.2e1 * pkin(3) * t48, -0.2e1 * pkin(3) * t46, 0.2e1 * t62 * qJ(4), t62 * qJ(4) ^ 2 + pkin(3) ^ 2, t26, 0.2e1 * t33 * t79, 0, t24, 0, 0, -t79 * t76, t33 * t76, -0.2e1 * t16 * t33 + 0.2e1 * t19 * t79, t16 ^ 2 + t19 ^ 2 + t40 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, t66, 0, t17, 0, 0, 0, 0, 0, 0, t9, t11, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, t46, 0, -pkin(3), 0, 0, 0, 0, 0, 0, -t79, t33, 0, t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, -t9, t31, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, -t33, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, t79, 0, t16, -t19, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t10;
