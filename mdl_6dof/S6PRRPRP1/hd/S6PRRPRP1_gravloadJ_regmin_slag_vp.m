% Calculate minimal parameter regressor of gravitation load for
% S6PRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% taug_reg [6x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:40
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRPRP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:38:48
% EndTime: 2021-01-16 02:38:50
% DurationCPUTime: 0.43s
% Computational Cost: add. (413->90), mult. (753->143), div. (0->0), fcn. (911->12), ass. (0->63)
t42 = sin(qJ(3));
t45 = cos(qJ(3));
t70 = cos(pkin(6));
t38 = sin(pkin(6));
t43 = sin(qJ(2));
t74 = t38 * t43;
t86 = -t42 * t74 + t70 * t45;
t37 = sin(pkin(10));
t62 = t37 * t70;
t69 = cos(pkin(10));
t82 = cos(qJ(2));
t19 = t43 * t62 - t69 * t82;
t73 = t38 * t45;
t85 = t19 * t42 + t37 * t73;
t55 = t70 * t69;
t21 = t37 * t82 + t43 * t55;
t36 = qJ(3) + pkin(11);
t34 = sin(t36);
t35 = cos(t36);
t61 = t38 * t69;
t10 = t21 * t34 + t35 * t61;
t75 = t37 * t38;
t12 = t19 * t34 + t35 * t75;
t16 = t34 * t74 - t70 * t35;
t51 = g(1) * t12 - g(2) * t10 - g(3) * t16;
t44 = cos(qJ(5));
t32 = t44 * pkin(5) + pkin(4);
t39 = -qJ(6) - pkin(9);
t57 = -t32 * t35 + t34 * t39;
t11 = t21 * t35 - t34 * t61;
t17 = t70 * t34 + t35 * t74;
t20 = t37 * t43 - t82 * t55;
t22 = t69 * t43 + t82 * t62;
t41 = sin(qJ(5));
t63 = t44 * t82;
t9 = t19 * t35 - t34 * t75;
t1 = -g(1) * (t22 * t44 + t41 * t9) - g(2) * (-t11 * t41 + t20 * t44) - g(3) * (-t17 * t41 - t38 * t63);
t83 = g(3) * t38;
t81 = t19 * t41;
t77 = t35 * t41;
t76 = t35 * t44;
t72 = t41 * t43;
t33 = t45 * pkin(3) + pkin(2);
t40 = qJ(4) + pkin(8);
t71 = -t19 * t40 - t22 * t33;
t65 = t38 * t82;
t66 = g(3) * (t33 * t65 + t40 * t74);
t64 = t41 * t82;
t59 = t85 * pkin(3);
t56 = t86 * pkin(3);
t52 = g(1) * t9 - g(2) * t11 - g(3) * t17;
t50 = -t21 * t42 - t45 * t61;
t49 = -g(1) * t19 + g(2) * t21 + g(3) * t74;
t47 = t50 * pkin(3);
t46 = g(1) * t22 + g(2) * t20 - g(3) * t65;
t14 = t20 * t33;
t7 = t46 * t34;
t6 = t51 * t44;
t5 = t51 * t41;
t4 = -g(1) * (-t22 * t76 - t81) - g(2) * (-t20 * t76 + t21 * t41) - (t35 * t63 + t72) * t83;
t3 = -g(1) * (-t19 * t44 + t22 * t77) - g(2) * (t20 * t77 + t21 * t44) - (-t35 * t64 + t43 * t44) * t83;
t2 = -g(1) * (-t22 * t41 + t44 * t9) - g(2) * (-t11 * t44 - t20 * t41) - g(3) * (-t17 * t44 + t38 * t64);
t8 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, t46, t49, 0, 0, 0, 0, 0, t46 * t45, -t46 * t42, t46 * t35, -t7, -t49, -g(1) * t71 - g(2) * (t21 * t40 - t14) - t66, 0, 0, 0, 0, 0, t4, t3, t4, t3, t7, -g(1) * (-pkin(5) * t81 + t57 * t22 + t71) - g(2) * (-t14 + (pkin(5) * t41 + t40) * t21 + t57 * t20) - t66 - (pkin(5) * t72 - t57 * t82) * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t85 - g(2) * t50 - g(3) * t86, -g(1) * (t19 * t45 - t42 * t75) - g(2) * (-t21 * t45 + t42 * t61) - g(3) * (-t70 * t42 - t43 * t73), -t51, -t52, 0, -g(1) * t59 - g(2) * t47 - g(3) * t56, 0, 0, 0, 0, 0, -t6, t5, -t6, t5, t52, -g(1) * (t12 * t32 + t9 * t39 + t59) - g(2) * (-t10 * t32 - t11 * t39 + t47) - g(3) * (-t16 * t32 - t17 * t39 + t56); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51;];
taug_reg = t8;
