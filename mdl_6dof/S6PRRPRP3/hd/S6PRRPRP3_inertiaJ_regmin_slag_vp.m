% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRPRP3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t45 = sin(pkin(11));
t47 = cos(pkin(11));
t49 = sin(qJ(5));
t73 = cos(qJ(5));
t81 = -t49 * t45 + t73 * t47;
t50 = sin(qJ(3));
t24 = t81 * t50;
t80 = -0.2e1 * t24;
t32 = t45 * t73 + t49 * t47;
t79 = -0.2e1 * t32;
t40 = -t47 * pkin(4) - pkin(3);
t78 = 0.2e1 * t40;
t52 = cos(qJ(3));
t77 = 0.2e1 * t52;
t76 = pkin(8) * t45;
t41 = t50 * pkin(8);
t75 = t52 * pkin(5);
t74 = t52 * pkin(8);
t34 = -t52 * pkin(3) - t50 * qJ(4) - pkin(2);
t29 = t47 * t34;
t66 = t47 * t50;
t12 = -pkin(9) * t66 + t29 + (-pkin(4) - t76) * t52;
t22 = t45 * t34 + t47 * t74;
t69 = t45 * t50;
t18 = -pkin(9) * t69 + t22;
t4 = t49 * t12 + t73 * t18;
t64 = pkin(9) + qJ(4);
t35 = t64 * t47;
t59 = t64 * t45;
t19 = t49 * t35 + t59 * t73;
t72 = t19 * t52;
t20 = t35 * t73 - t49 * t59;
t71 = t20 * t52;
t48 = cos(pkin(6));
t46 = sin(pkin(6));
t68 = t46 * sin(qJ(2));
t26 = -t48 * t52 + t50 * t68;
t70 = t26 * t32;
t67 = t46 * cos(qJ(2));
t33 = pkin(4) * t69 + t41;
t63 = t45 ^ 2 + t47 ^ 2;
t62 = t52 * qJ(6);
t23 = t32 * t50;
t27 = t48 * t50 + t52 * t68;
t16 = -t27 * t45 - t47 * t67;
t17 = t27 * t47 - t45 * t67;
t5 = -t16 * t73 + t49 * t17;
t60 = t26 * t23 + t5 * t52;
t3 = t73 * t12 - t49 * t18;
t6 = t49 * t16 + t17 * t73;
t58 = t26 * t24 + t6 * t52;
t57 = -pkin(3) * t50 + qJ(4) * t52;
t56 = -t16 * t45 + t17 * t47;
t21 = -t45 * t74 + t29;
t55 = -t21 * t45 + t22 * t47;
t44 = t50 ^ 2;
t25 = t26 ^ 2;
t15 = t26 * t81;
t10 = -pkin(5) * t81 - t32 * qJ(6) + t40;
t7 = t23 * pkin(5) - t24 * qJ(6) + t33;
t2 = -t3 + t75;
t1 = -t62 + t4;
t8 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16 ^ 2 + t17 ^ 2 + t25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 ^ 2 + t6 ^ 2 + t25; 0, 0, t67, -t68, 0, 0, 0, 0, 0, t52 * t67, -t50 * t67, -t16 * t52 + t26 * t69, t17 * t52 + t26 * t66 (-t16 * t47 - t17 * t45) * t50, t16 * t21 + t17 * t22 + t26 * t41, 0, 0, 0, 0, 0, t60, t58, t60, -t6 * t23 + t5 * t24, -t58, t6 * t1 + t5 * t2 + t26 * t7; 0, 1, 0, 0, t44, t50 * t77, 0, 0, 0, pkin(2) * t77, -0.2e1 * pkin(2) * t50, -0.2e1 * t21 * t52 + 0.2e1 * t44 * t76, 0.2e1 * t44 * pkin(8) * t47 + 0.2e1 * t22 * t52, 0.2e1 * (-t21 * t47 - t22 * t45) * t50, t44 * pkin(8) ^ 2 + t21 ^ 2 + t22 ^ 2, t24 ^ 2, t23 * t80, t52 * t80, t23 * t77, t52 ^ 2, 0.2e1 * t33 * t23 - 0.2e1 * t3 * t52, 0.2e1 * t33 * t24 + 0.2e1 * t4 * t52, 0.2e1 * t2 * t52 + 0.2e1 * t7 * t23, -0.2e1 * t1 * t23 + 0.2e1 * t2 * t24, -0.2e1 * t1 * t52 - 0.2e1 * t7 * t24, t1 ^ 2 + t2 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, -t27, -t26 * t47, t26 * t45, t56, -t26 * pkin(3) + qJ(4) * t56, 0, 0, 0, 0, 0, -t15, t70, -t15, t5 * t32 + t6 * t81, -t70, t26 * t10 + t5 * t19 + t6 * t20; 0, 0, 0, 0, 0, 0, t50, t52, 0, -t41, -t74, -pkin(8) * t66 + t45 * t57, pkin(8) * t69 + t47 * t57, t55, -pkin(3) * t41 + qJ(4) * t55, t24 * t32, -t32 * t23 + t24 * t81, -t32 * t52, -t81 * t52, 0, t40 * t23 - t33 * t81 + t72, t40 * t24 + t33 * t32 + t71, t10 * t23 - t7 * t81 + t72, t1 * t81 + t19 * t24 + t2 * t32 - t20 * t23, -t10 * t24 - t7 * t32 - t71, t1 * t20 + t7 * t10 + t2 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3) * t47, -0.2e1 * pkin(3) * t45, 0.2e1 * t63 * qJ(4), qJ(4) ^ 2 * t63 + pkin(3) ^ 2, t32 ^ 2, -t81 * t79, 0, 0, 0, -t81 * t78, t32 * t78, -0.2e1 * t10 * t81, 0.2e1 * t19 * t32 + 0.2e1 * t20 * t81, t10 * t79, t10 ^ 2 + t19 ^ 2 + t20 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, t66, 0, t41, 0, 0, 0, 0, 0, t23, t24, t23, 0, -t24, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, t45, 0, -pkin(3), 0, 0, 0, 0, 0, -t81, t32, -t81, 0, -t32, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t6, -t5, 0, t6, -t5 * pkin(5) + t6 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t23, -t52, t3, -t4, t3 - 0.2e1 * t75, -pkin(5) * t24 - t23 * qJ(6), -0.2e1 * t62 + t4, -t2 * pkin(5) + t1 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, t81, 0, -t19, -t20, -t19, -pkin(5) * t32 + qJ(6) * t81, t20, -t19 * pkin(5) + t20 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, t24, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t8;
