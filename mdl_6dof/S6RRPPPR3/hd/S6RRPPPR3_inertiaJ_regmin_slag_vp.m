% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPPPR3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t49 = cos(qJ(2));
t33 = t49 * qJ(4);
t65 = t49 * pkin(7);
t27 = -t33 + t65;
t71 = t27 ^ 2;
t44 = cos(pkin(9));
t45 = qJ(3) + pkin(4);
t30 = t44 * pkin(5) + t45;
t70 = -0.2e1 * t30;
t47 = sin(qJ(2));
t69 = -0.2e1 * t47;
t68 = 0.2e1 * t47;
t67 = -0.2e1 * t49;
t66 = 0.2e1 * t49;
t50 = -pkin(2) - pkin(3);
t38 = -qJ(5) + t50;
t64 = pkin(8) - t38;
t43 = sin(pkin(9));
t46 = sin(qJ(6));
t48 = cos(qJ(6));
t19 = t46 * t43 - t48 * t44;
t63 = t19 * t47;
t20 = t48 * t43 + t46 * t44;
t16 = t20 * t47;
t62 = t43 * t47;
t61 = t43 * t49;
t60 = t44 * t47;
t59 = t44 * t49;
t25 = -t49 * pkin(2) - t47 * qJ(3) - pkin(1);
t18 = t49 * pkin(3) - t25;
t12 = t47 * pkin(4) + t49 * qJ(5) + t18;
t35 = t47 * pkin(7);
t26 = -t47 * qJ(4) + t35;
t7 = t43 * t12 + t44 * t26;
t29 = t43 ^ 2 + t44 ^ 2;
t41 = t47 ^ 2;
t58 = t49 ^ 2 + t41;
t57 = t49 * qJ(3);
t6 = t44 * t12 - t43 * t26;
t56 = t7 * t43 + t6 * t44;
t3 = -t6 * t43 + t7 * t44;
t55 = -t47 * pkin(2) + t57;
t54 = -t38 * t47 - t45 * t49;
t52 = qJ(3) ^ 2;
t51 = 0.2e1 * qJ(3);
t24 = t64 * t44;
t23 = t64 * t43;
t17 = t29 * t38;
t15 = -t33 + (-pkin(5) * t43 + pkin(7)) * t49;
t14 = t19 * t49;
t13 = t20 * t49;
t9 = t46 * t23 - t48 * t24;
t8 = t48 * t23 + t46 * t24;
t5 = pkin(8) * t61 + t7;
t4 = t47 * pkin(5) + pkin(8) * t59 + t6;
t2 = t46 * t4 + t48 * t5;
t1 = t48 * t4 - t46 * t5;
t10 = [1, 0, 0, t41, t47 * t66, 0, 0, 0, pkin(1) * t66, pkin(1) * t69, t25 * t67, 0.2e1 * t58 * pkin(7), t25 * t69, t58 * pkin(7) ^ 2 + t25 ^ 2, t18 * t68, t18 * t67, -0.2e1 * t26 * t47 - 0.2e1 * t27 * t49, t18 ^ 2 + t26 ^ 2 + t71, -0.2e1 * t27 * t61 + 0.2e1 * t6 * t47, -0.2e1 * t27 * t59 - 0.2e1 * t7 * t47, t56 * t66, t6 ^ 2 + t7 ^ 2 + t71, t14 ^ 2, 0.2e1 * t14 * t13, t14 * t68, t13 * t68, t41, 0.2e1 * t1 * t47 - 0.2e1 * t15 * t13, 0.2e1 * t15 * t14 - 0.2e1 * t2 * t47; 0, 0, 0, 0, 0, t47, t49, 0, -t35, -t65, -t35, t55, t65, t55 * pkin(7), t27, t26, -t50 * t47 - t57, t27 * qJ(3) + t26 * t50, t27 * t44 + t54 * t43, -t27 * t43 + t54 * t44, -t3, t27 * t45 + t3 * t38, -t14 * t20, -t20 * t13 + t14 * t19, -t16, t63, 0, -t30 * t13 - t15 * t19 + t8 * t47, t30 * t14 - t15 * t20 - t9 * t47; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(2), 0, t51, pkin(2) ^ 2 + t52, t51, 0.2e1 * t50, 0, t50 ^ 2 + t52, 0.2e1 * t45 * t44, -0.2e1 * t45 * t43, -0.2e1 * t17, t29 * t38 ^ 2 + t45 ^ 2, t20 ^ 2, -0.2e1 * t20 * t19, 0, 0, 0, t19 * t70, t20 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, t35, 0, 0, -t47, t26, -t62, -t60, 0, t3, 0, 0, 0, 0, 0, -t16, t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(2), 0, 1, 0, t50, 0, 0, -t29, t17, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, t29, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, -t49, 0, t18, t60, -t62, t29 * t49, t56, 0, 0, 0, 0, 0, -t63, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t29, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, -t59, 0, t27, 0, 0, 0, 0, 0, -t13, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t43, 0, t45, 0, 0, 0, 0, 0, -t19, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t13, t47, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, t19, 0, t8, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t10;
