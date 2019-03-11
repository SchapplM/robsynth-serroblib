% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
% 
% Output:
% MM_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRRPR6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t52 = sin(qJ(3));
t55 = cos(qJ(4));
t40 = t55 * t52;
t50 = sin(qJ(6));
t51 = sin(qJ(4));
t54 = cos(qJ(6));
t73 = t54 * t51;
t20 = t50 * t40 - t52 * t73;
t90 = -0.2e1 * t20;
t69 = t51 * qJ(5);
t85 = pkin(4) + pkin(5);
t24 = t85 * t55 + pkin(3) + t69;
t89 = 0.2e1 * t24;
t61 = -t55 * pkin(4) - t69;
t32 = -pkin(3) + t61;
t88 = -0.2e1 * t32;
t87 = -0.2e1 * t52;
t56 = cos(qJ(3));
t86 = 0.2e1 * t56;
t84 = pkin(3) * t51;
t83 = pkin(3) * t55;
t82 = pkin(8) * t51;
t81 = pkin(8) * t55;
t66 = cos(pkin(6));
t49 = sin(pkin(6));
t78 = t49 * sin(qJ(2));
t22 = t52 * t78 - t66 * t56;
t80 = t22 * t51;
t79 = t22 * t55;
t77 = t49 * cos(qJ(2));
t76 = t51 * t52;
t75 = t51 * t55;
t74 = t51 * t56;
t72 = t55 * t56;
t33 = -t56 * pkin(3) - t52 * pkin(9) - pkin(2);
t71 = pkin(8) * t74 - t55 * t33;
t18 = pkin(8) * t72 + t51 * t33;
t45 = t51 ^ 2;
t47 = t55 ^ 2;
t70 = t45 + t47;
t68 = t55 * qJ(5);
t67 = t56 * qJ(5);
t65 = t52 * t86;
t44 = t56 * pkin(4);
t15 = t44 + t71;
t42 = t51 * pkin(9);
t64 = -t51 * pkin(10) + t42;
t23 = t66 * t52 + t56 * t78;
t9 = t23 * t51 + t55 * t77;
t63 = t22 * t76 + t9 * t56;
t14 = -t67 + t18;
t6 = t56 * pkin(5) - pkin(10) * t40 + t15;
t7 = pkin(10) * t76 + t14;
t1 = -t50 * t7 + t54 * t6;
t2 = t50 * t6 + t54 * t7;
t10 = t23 * t55 - t51 * t77;
t62 = t10 * t55 + t9 * t51;
t60 = -pkin(4) * t51 + t68;
t59 = t14 * t55 + t15 * t51;
t27 = t50 * t51 + t54 * t55;
t48 = t56 ^ 2;
t46 = t52 ^ 2;
t43 = t55 * pkin(9);
t37 = pkin(9) * t74;
t34 = -t55 * pkin(10) + t43;
t31 = t54 * qJ(5) - t50 * t85;
t30 = t50 * qJ(5) + t54 * t85;
t28 = -t50 * t55 + t73;
t21 = t27 * t52;
t19 = (pkin(8) - t60) * t52;
t13 = t54 * t34 + t50 * t64;
t12 = t50 * t34 - t54 * t64;
t11 = (-t85 * t51 - pkin(8) + t68) * t52;
t5 = t10 * t56 + t22 * t40;
t4 = t10 * t54 + t9 * t50;
t3 = t10 * t50 - t9 * t54;
t8 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 ^ 2 + t22 ^ 2 + t9 ^ 2, 0, 0, 0, 0, 0, 0, 0; 0, 0, t77, -t78, 0, 0, 0, 0, 0, t56 * t77, -t52 * t77, 0, 0, 0, 0, 0, t63, t5, t63 (-t10 * t51 + t55 * t9) * t52, -t5, t10 * t14 + t9 * t15 + t22 * t19, 0, 0, 0, 0, 0, -t22 * t20 - t3 * t56, -t22 * t21 - t4 * t56; 0, 1, 0, 0, t46, t65, 0, 0, 0, pkin(2) * t86, pkin(2) * t87, t47 * t46, -0.2e1 * t46 * t75, t72 * t87, t51 * t65, t48, 0.2e1 * t46 * t82 + 0.2e1 * t56 * t71, 0.2e1 * t18 * t56 + 0.2e1 * t46 * t81, 0.2e1 * t15 * t56 + 0.2e1 * t19 * t76, 0.2e1 * (-t14 * t51 + t15 * t55) * t52, -0.2e1 * t14 * t56 - 0.2e1 * t19 * t40, t14 ^ 2 + t15 ^ 2 + t19 ^ 2, t21 ^ 2, t21 * t90, t21 * t86, t56 * t90, t48, 0.2e1 * t1 * t56 + 0.2e1 * t11 * t20, 0.2e1 * t11 * t21 - 0.2e1 * t2 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, -t23, 0, 0, 0, 0, 0, -t79, t80, -t79, t62, -t80, t62 * pkin(9) + t22 * t32, 0, 0, 0, 0, 0, -t22 * t27, -t22 * t28; 0, 0, 0, 0, 0, 0, t52, t56, 0, -t52 * pkin(8), -t56 * pkin(8), t51 * t40 (-t45 + t47) * t52, -t74, -t72, 0, t37 + (-t81 - t84) * t52, pkin(9) * t72 + (t82 - t83) * t52, -t19 * t55 + t32 * t76 + t37, t59, -t19 * t51 + (-pkin(9) * t56 - t32 * t52) * t55, t59 * pkin(9) + t19 * t32, t21 * t28, -t28 * t20 - t21 * t27, t28 * t56, -t27 * t56, 0, t11 * t27 - t12 * t56 + t24 * t20, t11 * t28 - t13 * t56 + t24 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t45, 0.2e1 * t75, 0, 0, 0, 0.2e1 * t83, -0.2e1 * t84, t55 * t88, 0.2e1 * t70 * pkin(9), t51 * t88, t70 * pkin(9) ^ 2 + t32 ^ 2, t28 ^ 2, -0.2e1 * t28 * t27, 0, 0, 0, t27 * t89, t28 * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, -t10, -t9, 0, t10, -t9 * pkin(4) + t10 * qJ(5), 0, 0, 0, 0, 0, t3, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t76, -t56, -t71, -t18, -0.2e1 * t44 - t71, t61 * t52, -0.2e1 * t67 + t18, -t15 * pkin(4) + t14 * qJ(5), 0, 0, -t21, t20, -t56, -t30 * t56 - t1, -t31 * t56 + t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, t55, 0, -t42, -t43, -t42, t60, t43, t60 * pkin(9), 0, 0, -t28, t27, 0, t12, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t30, 0.2e1 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, t40, 0, t15, 0, 0, 0, 0, 0, t54 * t56, -t50 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, t42, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4), 0, 0, 0, 0, 0, -t54, t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t20, t56, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t27, 0, -t12, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t30, -t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, -t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t8;
