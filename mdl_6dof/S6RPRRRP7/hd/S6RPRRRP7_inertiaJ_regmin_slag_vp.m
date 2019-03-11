% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRRRP7
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
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRRP7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP7_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP7_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t55 = sin(pkin(10));
t56 = cos(pkin(10));
t59 = sin(qJ(3));
t82 = cos(qJ(3));
t34 = t59 * t55 - t82 * t56;
t88 = -0.2e1 * t34;
t87 = 0.2e1 * t34;
t57 = sin(qJ(5));
t60 = cos(qJ(4));
t58 = sin(qJ(4));
t81 = cos(qJ(5));
t67 = t81 * t58;
t38 = t57 * t60 + t67;
t86 = -0.2e1 * t38;
t45 = -t56 * pkin(2) - pkin(1);
t85 = 0.2e1 * t45;
t49 = -t60 * pkin(4) - pkin(3);
t84 = 0.2e1 * t49;
t83 = -pkin(9) - pkin(8);
t35 = t82 * t55 + t59 * t56;
t18 = t34 * pkin(3) - t35 * pkin(8) + t45;
t71 = pkin(7) + qJ(2);
t40 = t71 * t55;
t41 = t71 * t56;
t21 = -t59 * t40 + t82 * t41;
t73 = t60 * t21;
t10 = t73 + (-pkin(9) * t35 + t18) * t58;
t11 = t60 * t18 - t58 * t21;
t72 = t60 * t35;
t8 = t34 * pkin(4) - pkin(9) * t72 + t11;
t4 = t81 * t10 + t57 * t8;
t30 = t34 * pkin(5);
t50 = t57 * pkin(4);
t42 = t83 * t60;
t24 = -t57 * t42 - t83 * t67;
t80 = t24 * t34;
t77 = t57 * t58;
t25 = -t81 * t42 + t83 * t77;
t79 = t25 * t34;
t66 = t81 * t60;
t75 = t58 * t35;
t16 = t35 * t66 - t57 * t75;
t37 = -t66 + t77;
t78 = t37 * t16;
t22 = t37 * t34;
t23 = t38 * t34;
t76 = t58 * t34;
t74 = t58 * t60;
t70 = t55 ^ 2 + t56 ^ 2;
t69 = t35 * t88;
t29 = t34 * qJ(6);
t1 = t29 + t4;
t68 = t81 * pkin(4);
t3 = -t57 * t10 + t81 * t8;
t2 = -t3 - t30;
t65 = -pkin(3) * t35 - pkin(8) * t34;
t20 = t82 * t40 + t59 * t41;
t64 = -t37 * pkin(5) + t38 * qJ(6);
t14 = pkin(4) * t75 + t20;
t63 = 0.2e1 * pkin(5);
t61 = 0.2e1 * qJ(6);
t54 = t60 ^ 2;
t53 = t58 ^ 2;
t47 = t68 + pkin(5);
t44 = t50 + qJ(6);
t36 = t38 ^ 2;
t32 = t35 ^ 2;
t31 = t34 ^ 2;
t28 = t60 * t34;
t19 = t49 - t64;
t15 = t38 * t35;
t13 = t38 * t15;
t12 = t58 * t18 + t73;
t5 = t15 * pkin(5) - t16 * qJ(6) + t14;
t6 = [1, 0, 0, 0.2e1 * pkin(1) * t56, -0.2e1 * pkin(1) * t55, 0.2e1 * t70 * qJ(2), t70 * qJ(2) ^ 2 + pkin(1) ^ 2, t32, t69, 0, 0, 0, t34 * t85, t35 * t85, t54 * t32, -0.2e1 * t32 * t74, t72 * t87, t58 * t69, t31, 0.2e1 * t11 * t34 + 0.2e1 * t20 * t75, -0.2e1 * t12 * t34 + 0.2e1 * t20 * t72, t16 ^ 2, -0.2e1 * t16 * t15, t16 * t87, t15 * t88, t31, 0.2e1 * t14 * t15 + 0.2e1 * t3 * t34, 0.2e1 * t14 * t16 - 0.2e1 * t4 * t34, 0.2e1 * t5 * t15 - 0.2e1 * t2 * t34, -0.2e1 * t1 * t15 + 0.2e1 * t2 * t16, 0.2e1 * t1 * t34 - 0.2e1 * t5 * t16, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, -t56, t55, 0, -pkin(1), 0, 0, 0, 0, 0, t34, t35, 0, 0, 0, 0, 0, t28, -t76, 0, 0, 0, 0, 0, -t22, -t23, -t22, -t13 + t78, t23, t1 * t38 + t2 * t37; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37 ^ 2 + t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t34, 0, -t20, -t21, t58 * t72 (-t53 + t54) * t35, t76, t28, 0, -t20 * t60 + t65 * t58, t20 * t58 + t65 * t60, t16 * t38, -t13 - t78, t23, -t22, 0, t14 * t37 + t49 * t15 - t80, t14 * t38 + t49 * t16 - t79, t19 * t15 + t5 * t37 - t80, -t1 * t37 - t25 * t15 + t24 * t16 + t2 * t38, -t19 * t16 - t5 * t38 + t79, t1 * t25 + t5 * t19 + t2 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37 * t24 + t38 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t53, 0.2e1 * t74, 0, 0, 0, 0.2e1 * pkin(3) * t60, -0.2e1 * pkin(3) * t58, t36, t37 * t86, 0, 0, 0, t37 * t84, t38 * t84, 0.2e1 * t19 * t37, 0.2e1 * t24 * t38 - 0.2e1 * t25 * t37, t19 * t86, t19 ^ 2 + t24 ^ 2 + t25 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, -t75, t34, t11, -t12, 0, 0, t16, -t15, t34, t34 * t68 + t3, -t34 * t50 - t4, t47 * t34 - t2, -t44 * t15 - t47 * t16, t44 * t34 + t1, t1 * t44 - t2 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, -t58, 0, 0, 0, 0, 0, -t37, -t38, -t37, 0, t38, -t37 * t47 + t38 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, t60, 0, -t58 * pkin(8), -t60 * pkin(8), 0, 0, t38, -t37, 0, -t24, -t25, -t24, -t44 * t37 - t47 * t38, t25, -t24 * t47 + t25 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t68, -0.2e1 * t50, 0.2e1 * t47, 0, 0.2e1 * t44, t44 ^ 2 + t47 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t15, t34, t3, -t4, -t2 + t30, -pkin(5) * t16 - t15 * qJ(6), 0.2e1 * t29 + t4, -t2 * pkin(5) + t1 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, -t38, -t37, 0, t38, t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t37, 0, -t24, -t25, -t24, -pkin(5) * t38 - t37 * qJ(6), t25, -t24 * pkin(5) + t25 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t68, -t50, t63 + t68, 0, t61 + t50, t47 * pkin(5) + t44 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t63, 0, t61, pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, t16, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t6;
