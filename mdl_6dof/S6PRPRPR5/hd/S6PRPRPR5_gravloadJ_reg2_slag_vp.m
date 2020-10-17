% Calculate inertial parameters regressor of gravitation load for
% S6PRPRPR5
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
% Datum: 2019-03-08 19:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRPR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:59:06
% EndTime: 2019-05-04 22:59:07
% DurationCPUTime: 0.41s
% Computational Cost: add. (429->99), mult. (735->142), div. (0->0), fcn. (881->12), ass. (0->55)
t39 = pkin(11) + qJ(4);
t38 = cos(t39);
t37 = sin(t39);
t62 = qJ(5) * t37;
t42 = sin(pkin(6));
t48 = cos(qJ(2));
t67 = t42 * t48;
t75 = (pkin(4) * t38 + t62) * t67;
t74 = g(3) * t42;
t41 = sin(pkin(10));
t46 = sin(qJ(2));
t60 = cos(pkin(10));
t61 = cos(pkin(6));
t51 = t61 * t60;
t24 = t41 * t46 - t48 * t51;
t73 = t24 * t38;
t57 = t41 * t61;
t26 = t60 * t46 + t48 * t57;
t72 = t26 * t38;
t45 = sin(qJ(6));
t71 = t37 * t45;
t47 = cos(qJ(6));
t70 = t37 * t47;
t69 = t41 * t42;
t68 = t42 * t46;
t66 = t45 * t48;
t65 = t47 * t48;
t25 = t41 * t48 + t46 * t51;
t43 = cos(pkin(11));
t36 = t43 * pkin(3) + pkin(2);
t44 = -pkin(8) - qJ(3);
t64 = -t24 * t36 - t25 * t44;
t27 = -t46 * t57 + t60 * t48;
t63 = -t26 * t36 - t27 * t44;
t56 = t42 * t60;
t10 = t25 * t38 - t37 * t56;
t9 = t25 * t37 + t38 * t56;
t59 = -t9 * pkin(4) + t10 * qJ(5);
t11 = t27 * t37 - t38 * t69;
t12 = t27 * t38 + t37 * t69;
t58 = -t11 * pkin(4) + t12 * qJ(5);
t20 = t37 * t68 - t61 * t38;
t21 = t61 * t37 + t38 * t68;
t55 = -t20 * pkin(4) + t21 * qJ(5);
t54 = -pkin(4) * t73 - t24 * t62 + t64;
t53 = -pkin(4) * t72 - t26 * t62 + t63;
t29 = t36 * t67;
t52 = -t44 * t68 + t29;
t2 = g(1) * t11 + g(2) * t9 + g(3) * t20;
t50 = g(1) * t12 + g(2) * t10 + g(3) * t21;
t5 = -g(1) * t26 - g(2) * t24 + g(3) * t67;
t49 = g(1) * t27 + g(2) * t25 + g(3) * t68;
t4 = t5 * t38;
t3 = t5 * t37;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, t49, 0, 0, 0, 0, 0, 0, 0, 0, -t5 * t43, t5 * sin(pkin(11)) -t49, -g(1) * (-t26 * pkin(2) + t27 * qJ(3)) - g(2) * (-t24 * pkin(2) + t25 * qJ(3)) - (pkin(2) * t48 + qJ(3) * t46) * t74, 0, 0, 0, 0, 0, 0, -t4, t3, -t49, -g(1) * t63 - g(2) * t64 - g(3) * t52, 0, 0, 0, 0, 0, 0, -t49, t4, -t3, -g(1) * t53 - g(2) * t54 - g(3) * (t52 + t75) 0, 0, 0, 0, 0, 0, -g(1) * (-t26 * t71 + t27 * t47) - g(2) * (-t24 * t71 + t25 * t47) - (t37 * t66 + t46 * t47) * t74, -g(1) * (-t26 * t70 - t27 * t45) - g(2) * (-t24 * t70 - t25 * t45) - (t37 * t65 - t45 * t46) * t74, -t4, -g(1) * (t27 * pkin(5) - pkin(9) * t72 + t53) - g(2) * (t25 * pkin(5) - pkin(9) * t73 + t54) - g(3) * (t29 + t75) - (pkin(9) * t38 * t48 + (pkin(5) - t44) * t46) * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t50, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, -t50, -g(1) * t58 - g(2) * t59 - g(3) * t55, 0, 0, 0, 0, 0, 0, -t50 * t45, -t50 * t47, t2, -g(1) * (-t11 * pkin(9) + t58) - g(2) * (-t9 * pkin(9) + t59) - g(3) * (-t20 * pkin(9) + t55); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t11 * t47 - t26 * t45) - g(2) * (-t24 * t45 + t9 * t47) - g(3) * (t20 * t47 + t42 * t66) -g(1) * (-t11 * t45 - t26 * t47) - g(2) * (-t24 * t47 - t9 * t45) - g(3) * (-t20 * t45 + t42 * t65) 0, 0;];
taug_reg  = t1;
