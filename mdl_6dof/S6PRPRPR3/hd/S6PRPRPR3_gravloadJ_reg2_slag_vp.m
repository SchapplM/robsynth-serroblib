% Calculate inertial parameters regressor of gravitation load for
% S6PRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRPR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t46 = sin(pkin(11));
t53 = sin(qJ(2));
t56 = cos(qJ(2));
t74 = cos(pkin(11));
t36 = -t56 * t46 - t53 * t74;
t47 = sin(pkin(10));
t49 = cos(pkin(10));
t50 = cos(pkin(6));
t64 = -t53 * t46 + t56 * t74;
t60 = t64 * t50;
t19 = t47 * t36 + t49 * t60;
t52 = sin(qJ(4));
t75 = qJ(5) * t52;
t55 = cos(qJ(4));
t92 = t19 * t55;
t96 = pkin(4) * t92 + t19 * t75;
t22 = t49 * t36 - t47 * t60;
t91 = t22 * t55;
t95 = pkin(4) * t91 + t22 * t75;
t48 = sin(pkin(6));
t33 = t64 * t48;
t90 = t33 * t55;
t94 = pkin(4) * t90 + t33 * t75;
t77 = t36 * t50;
t23 = t47 * t77 + t49 * t64;
t18 = -t47 * t64 + t49 * t77;
t93 = pkin(5) + pkin(8);
t88 = t47 * t53;
t87 = t48 * t52;
t86 = t48 * t55;
t85 = t48 * t56;
t83 = t50 * t53;
t82 = t50 * t56;
t51 = sin(qJ(6));
t81 = t51 * t52;
t54 = cos(qJ(6));
t80 = t52 * t54;
t76 = pkin(2) * t85 + t33 * pkin(3);
t73 = t49 * t82;
t8 = -t18 * t52 + t49 * t86;
t9 = -t18 * t55 - t49 * t87;
t72 = -t8 * pkin(4) + t9 * qJ(5);
t10 = t23 * t52 - t47 * t86;
t11 = t23 * t55 + t47 * t87;
t71 = -t10 * pkin(4) + t11 * qJ(5);
t34 = t36 * t48;
t25 = -t34 * t52 - t50 * t55;
t26 = -t34 * t55 + t50 * t52;
t68 = -t25 * pkin(4) + t26 * qJ(5);
t67 = -t34 * pkin(8) + t76;
t37 = pkin(2) * t73;
t66 = -pkin(2) * t88 + t19 * pkin(3) + t37;
t65 = -t47 * t82 - t49 * t53;
t2 = g(1) * t10 + g(2) * t8 + g(3) * t25;
t63 = g(1) * t11 + g(2) * t9 + g(3) * t26;
t5 = -g(1) * t23 + g(2) * t18 + g(3) * t34;
t62 = g(1) * t22 + g(2) * t19 + g(3) * t33;
t61 = -t18 * pkin(8) + t66;
t59 = pkin(2) * t65 + t22 * pkin(3);
t58 = -g(1) * t65 - g(3) * t85;
t57 = pkin(8) * t23 + t59;
t32 = -g(3) * t50 + (-g(1) * t47 + g(2) * t49) * t48;
t4 = t62 * t55;
t3 = t62 * t52;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2) * (t73 - t88) + t58, -g(1) * (t47 * t83 - t49 * t56) - g(2) * (-t47 * t56 - t49 * t83) + g(3) * t48 * t53, 0, 0, 0, 0, 0, 0, 0, 0, -t62, -t5, 0, -g(2) * t37 + (g(2) * t88 + t58) * pkin(2), 0, 0, 0, 0, 0, 0, -t4, t3, t5, -g(1) * t57 - g(2) * t61 - g(3) * t67, 0, 0, 0, 0, 0, 0, t5, t4, -t3, -g(1) * (t57 + t95) - g(2) * (t61 + t96) - g(3) * (t67 + t94) 0, 0, 0, 0, 0, 0, -g(1) * (t22 * t81 + t23 * t54) - g(2) * (-t18 * t54 + t19 * t81) - g(3) * (t33 * t81 - t34 * t54) -g(1) * (t22 * t80 - t23 * t51) - g(2) * (t18 * t51 + t19 * t80) - g(3) * (t33 * t80 + t34 * t51) -t4, -g(1) * (pkin(9) * t91 + t23 * t93 + t59 + t95) - g(2) * (pkin(9) * t92 - t18 * t93 + t66 + t96) - g(3) * (pkin(9) * t90 - t34 * t93 + t76 + t94); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t63, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, -t63, -g(1) * t71 - g(2) * t72 - g(3) * t68, 0, 0, 0, 0, 0, 0, -t63 * t51, -t63 * t54, t2, -g(1) * (-t10 * pkin(9) + t71) - g(2) * (-t8 * pkin(9) + t72) - g(3) * (-t25 * pkin(9) + t68); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t10 * t54 + t22 * t51) - g(2) * (t19 * t51 + t8 * t54) - g(3) * (t25 * t54 + t33 * t51) -g(1) * (-t10 * t51 + t22 * t54) - g(2) * (t19 * t54 - t8 * t51) - g(3) * (-t25 * t51 + t33 * t54) 0, 0;];
taug_reg  = t1;
