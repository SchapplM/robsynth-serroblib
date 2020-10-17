% Calculate inertial parameters regressor of gravitation load for
% S6PRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRPPR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 02:46:34
% EndTime: 2019-05-05 02:46:35
% DurationCPUTime: 0.52s
% Computational Cost: add. (462->114), mult. (827->167), div. (0->0), fcn. (984->12), ass. (0->63)
t46 = sin(qJ(3));
t49 = cos(qJ(3));
t70 = cos(pkin(6));
t43 = sin(pkin(6));
t47 = sin(qJ(2));
t78 = t43 * t47;
t88 = -t46 * t78 + t49 * t70;
t50 = cos(qJ(2));
t42 = sin(pkin(10));
t66 = t42 * t70;
t69 = cos(pkin(10));
t27 = -t47 * t66 + t50 * t69;
t77 = t43 * t49;
t87 = -t27 * t46 + t42 * t77;
t41 = qJ(3) + pkin(11);
t40 = cos(t41);
t39 = sin(t41);
t71 = qJ(5) * t39;
t76 = t43 * t50;
t86 = (pkin(4) * t40 + t71) * t76;
t85 = g(3) * t43;
t58 = t70 * t69;
t24 = t42 * t47 - t50 * t58;
t84 = t24 * t40;
t26 = t47 * t69 + t50 * t66;
t83 = t26 * t40;
t45 = sin(qJ(6));
t81 = t39 * t45;
t48 = cos(qJ(6));
t80 = t39 * t48;
t79 = t42 * t43;
t75 = t45 * t50;
t74 = t48 * t50;
t25 = t42 * t50 + t47 * t58;
t38 = pkin(3) * t49 + pkin(2);
t44 = -qJ(4) - pkin(8);
t73 = -t24 * t38 - t25 * t44;
t72 = -t26 * t38 - t27 * t44;
t65 = t43 * t69;
t63 = -pkin(4) * t84 - t24 * t71 + t73;
t62 = -pkin(4) * t83 - t26 * t71 + t72;
t61 = t87 * pkin(3);
t30 = t38 * t76;
t60 = -t44 * t78 + t30;
t59 = t88 * pkin(3);
t11 = t27 * t39 - t40 * t79;
t20 = t39 * t78 - t40 * t70;
t9 = t25 * t39 + t40 * t65;
t2 = g(1) * t11 + g(2) * t9 + g(3) * t20;
t10 = t25 * t40 - t39 * t65;
t12 = t27 * t40 + t39 * t79;
t21 = t39 * t70 + t40 * t78;
t57 = g(1) * t12 + g(2) * t10 + g(3) * t21;
t56 = -t25 * t46 - t49 * t65;
t55 = -pkin(4) * t11 + qJ(5) * t12 + t61;
t5 = -g(1) * t26 - g(2) * t24 + g(3) * t76;
t54 = g(1) * t27 + g(2) * t25 + g(3) * t78;
t53 = -pkin(4) * t20 + qJ(5) * t21 + t59;
t52 = t56 * pkin(3);
t51 = -pkin(4) * t9 + t10 * qJ(5) + t52;
t4 = t5 * t40;
t3 = t5 * t39;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, t54, 0, 0, 0, 0, 0, 0, 0, 0, -t5 * t49, t5 * t46, -t54, -g(1) * (-pkin(2) * t26 + pkin(8) * t27) - g(2) * (-pkin(2) * t24 + pkin(8) * t25) - (pkin(2) * t50 + pkin(8) * t47) * t85, 0, 0, 0, 0, 0, 0, -t4, t3, -t54, -g(1) * t72 - g(2) * t73 - g(3) * t60, 0, 0, 0, 0, 0, 0, -t54, t4, -t3, -g(1) * t62 - g(2) * t63 - g(3) * (t60 + t86) 0, 0, 0, 0, 0, 0, -g(1) * (-t26 * t81 + t27 * t48) - g(2) * (-t24 * t81 + t25 * t48) - (t39 * t75 + t47 * t48) * t85, -g(1) * (-t26 * t80 - t27 * t45) - g(2) * (-t24 * t80 - t25 * t45) - (t39 * t74 - t45 * t47) * t85, -t4, -g(1) * (pkin(5) * t27 - pkin(9) * t83 + t62) - g(2) * (pkin(5) * t25 - pkin(9) * t84 + t63) - g(3) * (t30 + t86) - (pkin(9) * t40 * t50 + (pkin(5) - t44) * t47) * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t87 - g(2) * t56 - g(3) * t88, -g(1) * (-t27 * t49 - t46 * t79) - g(2) * (-t25 * t49 + t46 * t65) - g(3) * (-t46 * t70 - t47 * t77) 0, 0, 0, 0, 0, 0, 0, 0, t2, t57, 0, -g(1) * t61 - g(2) * t52 - g(3) * t59, 0, 0, 0, 0, 0, 0, 0, -t2, -t57, -g(1) * t55 - g(2) * t51 - g(3) * t53, 0, 0, 0, 0, 0, 0, -t57 * t45, -t57 * t48, t2, -g(1) * (-pkin(9) * t11 + t55) - g(2) * (-t9 * pkin(9) + t51) - g(3) * (-pkin(9) * t20 + t53); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t11 * t48 - t26 * t45) - g(2) * (-t24 * t45 + t48 * t9) - g(3) * (t20 * t48 + t43 * t75) -g(1) * (-t11 * t45 - t26 * t48) - g(2) * (-t24 * t48 - t45 * t9) - g(3) * (-t20 * t45 + t43 * t74) 0, 0;];
taug_reg  = t1;
