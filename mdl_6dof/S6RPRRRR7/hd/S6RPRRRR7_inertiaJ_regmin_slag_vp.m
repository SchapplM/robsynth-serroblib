% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% 
% Output:
% MM_reg [((6+1)*6/2)x34]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRRR7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR7_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR7_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t49 = sin(qJ(4));
t50 = sin(qJ(3));
t52 = cos(qJ(4));
t53 = cos(qJ(3));
t29 = t49 * t53 + t52 * t50;
t30 = -t49 * t50 + t52 * t53;
t48 = sin(qJ(5));
t70 = cos(qJ(5));
t15 = t48 * t29 - t70 * t30;
t14 = t15 ^ 2;
t47 = sin(qJ(6));
t80 = t15 * t47;
t51 = cos(qJ(6));
t79 = t15 * t51;
t56 = t70 * t29 + t48 * t30;
t78 = t56 ^ 2;
t11 = t47 * t56;
t12 = t51 * t56;
t36 = t50 * pkin(3) + qJ(2);
t22 = t29 * pkin(4) + t36;
t77 = 0.2e1 * t22;
t76 = 0.2e1 * t36;
t75 = 0.2e1 * qJ(2);
t74 = pkin(5) * t47;
t54 = -pkin(1) - pkin(7);
t31 = (-pkin(8) + t54) * t50;
t40 = t53 * t54;
t32 = -t53 * pkin(8) + t40;
t19 = -t49 * t31 + t52 * t32;
t55 = -t30 * pkin(9) + t19;
t20 = -t52 * t31 - t49 * t32;
t9 = -t29 * pkin(9) - t20;
t4 = t48 * t9 - t70 * t55;
t73 = t4 * t51;
t72 = t48 * pkin(4);
t71 = t49 * pkin(3);
t43 = t52 * pkin(3);
t39 = t43 + pkin(4);
t61 = -t70 * t39 + t48 * t71;
t23 = -pkin(5) + t61;
t68 = t23 * t51;
t42 = t70 * pkin(4);
t38 = -t42 - pkin(5);
t67 = t38 * t51;
t65 = t47 * t51;
t63 = 0.2e1 * t15 * t56;
t62 = t70 * t71;
t60 = pkin(5) * t15 - pkin(10) * t56;
t59 = -t78 - t14;
t26 = -t48 * t39 - t62;
t24 = pkin(10) - t26;
t58 = -t15 * t23 - t24 * t56;
t37 = pkin(10) + t72;
t57 = -t15 * t38 - t37 * t56;
t46 = t51 ^ 2;
t45 = t47 ^ 2;
t44 = pkin(5) * t51;
t35 = 0.2e1 * t65;
t34 = t38 * t47;
t21 = t23 * t47;
t10 = t47 * t79;
t7 = (-t45 + t46) * t15;
t6 = pkin(5) * t56 + pkin(10) * t15 + t22;
t5 = t48 * t55 + t70 * t9;
t3 = t4 * t47;
t2 = t47 * t6 + t51 * t5;
t1 = -t47 * t5 + t51 * t6;
t8 = [1, 0, 0, -2 * pkin(1), t75 (pkin(1) ^ 2) + qJ(2) ^ 2, t53 ^ 2, -0.2e1 * t53 * t50, 0, 0, 0, t50 * t75, t53 * t75, t30 ^ 2, -0.2e1 * t30 * t29, 0, 0, 0, t29 * t76, t30 * t76, t14, t63, 0, 0, 0, t56 * t77, -t15 * t77, t46 * t14, -0.2e1 * t14 * t65, -0.2e1 * t56 * t79, t47 * t63, t78, 0.2e1 * t1 * t56 - 0.2e1 * t4 * t80, -0.2e1 * t2 * t56 - 0.2e1 * t4 * t79; 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59 * t47, t59 * t51; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t53, -t50, 0, t40, -t50 * t54, 0, 0, t30, -t29, 0, t19, t20, 0, 0, -t15, -t56, 0, -t4, -t5, -t10, -t7, t11, t12, 0, t58 * t47 - t73, t58 * t51 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, -t50, 0, 0, 0, 0, 0, t30, -t29, 0, 0, 0, 0, 0, -t15, -t56, 0, 0, 0, 0, 0, -t79, t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t43, -0.2e1 * t71, 0, 0, 0, 0, 1, -0.2e1 * t61, 0.2e1 * t26, t45, t35, 0, 0, 0, -0.2e1 * t68, 0.2e1 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t29, 0, t19, t20, 0, 0, -t15, -t56, 0, -t4, -t5, -t10, -t7, t11, t12, 0, t57 * t47 - t73, t57 * t51 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t29, 0, 0, 0, 0, 0, -t15, -t56, 0, 0, 0, 0, 0, -t79, t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t43, -t71, 0, 0, 0, 0, 1, t42 - t61, -t62 + (-pkin(4) - t39) * t48, t45, t35, 0, 0, 0 (-t23 - t38) * t51, t34 + t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t42, -0.2e1 * t72, t45, t35, 0, 0, 0, -0.2e1 * t67, 0.2e1 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t56, 0, -t4, -t5, -t10, -t7, t11, t12, 0, t60 * t47 - t73, t60 * t51 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t56, 0, 0, 0, 0, 0, -t79, t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -t61, t26, t45, t35, 0, 0, 0, t44 - t68, t21 - t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t42, -t72, t45, t35, 0, 0, 0, t44 - t67, t34 - t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t45, t35, 0, 0, 0, 0.2e1 * t44, -0.2e1 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, t80, t56, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, t51, 0, -t47 * t24, -t51 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, t51, 0, -t47 * t37, -t51 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, t51, 0, -t47 * pkin(10), -t51 * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t8;
