% Calculate inertial parameters regressor of gravitation load for
% S6PRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRPPR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 03:27:40
% EndTime: 2019-05-05 03:27:42
% DurationCPUTime: 0.62s
% Computational Cost: add. (410->126), mult. (975->169), div. (0->0), fcn. (1180->12), ass. (0->67)
t43 = sin(qJ(3));
t45 = cos(qJ(3));
t84 = pkin(3) * t45 + qJ(4) * t43;
t44 = sin(qJ(2));
t46 = cos(qJ(2));
t66 = cos(pkin(10));
t67 = cos(pkin(6));
t52 = t67 * t66;
t65 = sin(pkin(10));
t20 = t44 * t65 - t46 * t52;
t83 = t84 * t20;
t51 = t67 * t65;
t22 = t44 * t66 + t46 * t51;
t82 = t84 * t22;
t81 = pkin(4) + pkin(8);
t40 = sin(pkin(6));
t79 = g(3) * t40;
t41 = cos(pkin(11));
t35 = t41 * pkin(5) + pkin(4);
t78 = pkin(8) + t35;
t38 = pkin(11) + qJ(6);
t36 = sin(t38);
t77 = t36 * t43;
t37 = cos(t38);
t76 = t37 * t43;
t39 = sin(pkin(11));
t75 = t39 * t43;
t74 = t40 * t44;
t73 = t40 * t46;
t72 = t41 * t43;
t71 = t43 * t46;
t70 = pkin(2) * t73 + pkin(8) * t74;
t68 = qJ(5) * t45;
t17 = t20 * pkin(2);
t64 = -t17 - t83;
t18 = t22 * pkin(2);
t63 = -t18 - t82;
t21 = t44 * t52 + t46 * t65;
t62 = t21 * pkin(8) - t17;
t23 = -t44 * t51 + t46 * t66;
t61 = t23 * pkin(8) - t18;
t60 = pkin(5) * t39 + qJ(4);
t57 = t40 * t66;
t10 = t21 * t45 - t43 * t57;
t9 = t21 * t43 + t45 * t57;
t7 = t9 * pkin(3);
t59 = t10 * qJ(4) - t7;
t56 = t40 * t65;
t12 = t23 * t45 + t43 * t56;
t11 = t23 * t43 - t45 * t56;
t8 = t11 * pkin(3);
t58 = t12 * qJ(4) - t8;
t24 = t43 * t74 - t45 * t67;
t19 = t24 * pkin(3);
t25 = t43 * t67 + t45 * t74;
t55 = t25 * qJ(4) - t19;
t54 = t84 * t73 + t70;
t53 = g(3) * t54;
t42 = -pkin(9) - qJ(5);
t50 = pkin(5) * t75 - t42 * t45;
t2 = g(1) * t11 + g(2) * t9 + g(3) * t24;
t49 = g(1) * t12 + g(2) * t10 + g(3) * t25;
t48 = -g(1) * t22 - g(2) * t20 + g(3) * t73;
t47 = g(1) * t23 + g(2) * t21 + g(3) * t74;
t5 = t48 * t45;
t4 = t48 * t43;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, t47, 0, 0, 0, 0, 0, 0, 0, 0, -t5, t4, -t47, -g(1) * t61 - g(2) * t62 - g(3) * t70, 0, 0, 0, 0, 0, 0, -t47, t5, -t4, -g(1) * (t61 - t82) - g(2) * (t62 - t83) - t53, 0, 0, 0, 0, 0, 0, -g(1) * (-t22 * t75 + t23 * t41) - g(2) * (-t20 * t75 + t21 * t41) - (t39 * t71 + t41 * t44) * t79, -g(1) * (-t22 * t72 - t23 * t39) - g(2) * (-t20 * t72 - t21 * t39) - (-t39 * t44 + t41 * t71) * t79, -t5, -g(1) * (-t22 * t68 + t23 * t81 + t63) - g(2) * (-t20 * t68 + t21 * t81 + t64) - g(3) * ((pkin(4) * t44 + t46 * t68) * t40 + t54) 0, 0, 0, 0, 0, 0, -g(1) * (-t22 * t77 + t23 * t37) - g(2) * (-t20 * t77 + t21 * t37) - (t36 * t71 + t37 * t44) * t79, -g(1) * (-t22 * t76 - t23 * t36) - g(2) * (-t20 * t76 - t21 * t36) - (-t36 * t44 + t37 * t71) * t79, -t5, -g(1) * (-t22 * t50 + t23 * t78 + t63) - g(2) * (-t20 * t50 + t21 * t78 + t64) - t53 - (t35 * t44 + t46 * t50) * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t49, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, -t49, -g(1) * t58 - g(2) * t59 - g(3) * t55, 0, 0, 0, 0, 0, 0, -t49 * t39, -t49 * t41, t2, -g(1) * (-t11 * qJ(5) + t58) - g(2) * (-t9 * qJ(5) + t59) - g(3) * (-t24 * qJ(5) + t55) 0, 0, 0, 0, 0, 0, -t49 * t36, -t49 * t37, t2, -g(1) * (t11 * t42 + t12 * t60 - t8) - g(2) * (t10 * t60 + t9 * t42 - t7) - g(3) * (t24 * t42 + t25 * t60 - t19); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t11 * t37 - t22 * t36) - g(2) * (-t20 * t36 + t9 * t37) - g(3) * (t24 * t37 + t36 * t73) -g(1) * (-t11 * t36 - t22 * t37) - g(2) * (-t20 * t37 - t9 * t36) - g(3) * (-t24 * t36 + t37 * t73) 0, 0;];
taug_reg  = t1;
