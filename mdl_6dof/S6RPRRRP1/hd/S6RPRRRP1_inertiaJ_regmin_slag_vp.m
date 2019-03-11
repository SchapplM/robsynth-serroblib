% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% MM_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRRP1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t49 = sin(qJ(4));
t83 = t49 * pkin(3);
t40 = pkin(9) + t83;
t48 = sin(qJ(5));
t44 = t48 ^ 2;
t51 = cos(qJ(5));
t45 = t51 ^ 2;
t66 = t44 + t45;
t93 = t66 * t40;
t50 = sin(qJ(3));
t77 = cos(qJ(4));
t78 = cos(qJ(3));
t29 = t49 * t78 + t50 * t77;
t56 = pkin(5) * t48 - t51 * qJ(6);
t92 = t56 * t29;
t47 = cos(pkin(10));
t38 = -t47 * pkin(1) - pkin(2);
t30 = -pkin(3) * t78 + t38;
t91 = 0.2e1 * t30;
t90 = -0.2e1 * t48;
t89 = 0.2e1 * t50;
t88 = -0.2e1 * t51;
t28 = t49 * t50 - t77 * t78;
t85 = t29 * pkin(9);
t11 = t28 * pkin(4) + t30 - t85;
t46 = sin(pkin(10));
t37 = t46 * pkin(1) + pkin(7);
t24 = (-pkin(8) - t37) * t50;
t62 = t78 * t37;
t25 = pkin(8) * t78 + t62;
t13 = t49 * t24 + t25 * t77;
t5 = t48 * t11 + t51 * t13;
t87 = pkin(9) * t28;
t86 = t28 * pkin(5);
t84 = t48 * pkin(9);
t82 = t51 * pkin(9);
t12 = -t77 * t24 + t49 * t25;
t6 = t12 + t92;
t81 = t6 * t48;
t80 = t6 * t51;
t63 = t77 * pkin(3);
t41 = -t63 - pkin(4);
t79 = pkin(4) - t41;
t76 = t12 * t51;
t75 = t28 * t40;
t74 = t44 * t29;
t73 = t48 * t29;
t72 = t48 * t40;
t71 = t48 * t51;
t22 = t51 * t29;
t70 = t51 * t40;
t57 = -t51 * pkin(5) - t48 * qJ(6);
t31 = -pkin(4) + t57;
t23 = -t63 + t31;
t69 = -t23 - t31;
t67 = pkin(9) * t66;
t65 = t28 * qJ(6);
t64 = -0.2e1 * t29 * t28;
t61 = -t51 * t11 + t48 * t13;
t59 = -pkin(4) * t29 - t87;
t2 = t65 + t5;
t3 = t61 - t86;
t1 = t2 * t51 + t3 * t48;
t58 = -t29 * t31 + t87;
t55 = t23 * t29 - t75;
t54 = t29 * t41 - t75;
t35 = 0.2e1 * t71;
t27 = t29 ^ 2;
t26 = t28 ^ 2;
t21 = t51 * t28;
t20 = t45 * t29;
t19 = t45 * t27;
t18 = t48 * t28;
t16 = t48 * t22;
t15 = t20 + t74;
t14 = t20 - t74;
t10 = t12 * t48;
t4 = [1, 0, 0 (t46 ^ 2 + t47 ^ 2) * pkin(1) ^ 2, t50 ^ 2, t78 * t89, 0, 0, 0, -0.2e1 * t38 * t78, t38 * t89, t27, t64, 0, 0, 0, t28 * t91, t29 * t91, t19, -0.2e1 * t27 * t71, 0.2e1 * t28 * t22, t48 * t64, t26, 0.2e1 * t12 * t73 - 0.2e1 * t28 * t61, 0.2e1 * t12 * t22 - 0.2e1 * t5 * t28, -0.2e1 * t3 * t28 + 0.2e1 * t6 * t73, 0.2e1 * (-t2 * t48 + t3 * t51) * t29, 0.2e1 * t2 * t28 - 0.2e1 * t22 * t6, t2 ^ 2 + t3 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 * t29 + t6 * t28; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44 * t27 + t19 + t26; 0, 0, 0, 0, 0, 0, t50, t78, 0, -t50 * t37, -t62, 0, 0, t29, -t28, 0, -t12, -t13, t16, t14, t18, t21, 0, t48 * t54 - t76, t51 * t54 + t10, t48 * t55 - t80, t1, -t51 * t55 - t81, t1 * t40 + t6 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, -t50, 0, 0, 0, 0, 0, -t28, -t29, 0, 0, 0, 0, 0, -t21, t18, -t21, t15, -t18, t28 * t23 + t29 * t93; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t63, -0.2e1 * t83, t44, t35, 0, 0, 0, t41 * t88, 0.2e1 * t41 * t48, t23 * t88, 0.2e1 * t93, t23 * t90, t40 ^ 2 * t66 + t23 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t28, 0, -t12, -t13, t16, t14, t18, t21, 0, t48 * t59 - t76, t51 * t59 + t10, -t48 * t58 - t80, t1, t51 * t58 - t81, pkin(9) * t1 + t6 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, -t29, 0, 0, 0, 0, 0, -t21, t18, -t21, t15, -t18, t28 * t31 + t66 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t63, -t83, t44, t35, 0, 0, 0, t79 * t51, -t79 * t48, t69 * t51, t67 + t93, t69 * t48, pkin(9) * t93 + t23 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t44, t35, 0, 0, 0, 0.2e1 * pkin(4) * t51, pkin(4) * t90, t31 * t88, 0.2e1 * t67, t31 * t90, pkin(9) ^ 2 * t66 + t31 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t73, t28, -t61, -t5, -t61 + 0.2e1 * t86, t57 * t29, 0.2e1 * t65 + t5, -t3 * pkin(5) + t2 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73, -t22, -t73, 0, t22, -t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t51, 0, -t72, -t70, -t72, -t56, t70, -t56 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t51, 0, -t84, -t82, -t84, -t56, t82, -t56 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, t22, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, 0, t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, 0, t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t4;
