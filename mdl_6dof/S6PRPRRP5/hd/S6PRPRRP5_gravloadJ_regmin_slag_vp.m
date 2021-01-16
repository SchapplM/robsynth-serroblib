% Calculate minimal parameter regressor of gravitation load for
% S6PRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% taug_reg [6x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 01:52
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRRP5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 01:51:50
% EndTime: 2021-01-16 01:51:52
% DurationCPUTime: 0.45s
% Computational Cost: add. (256->106), mult. (662->169), div. (0->0), fcn. (799->10), ass. (0->75)
t38 = sin(pkin(6));
t88 = g(3) * t38;
t37 = sin(pkin(10));
t39 = cos(pkin(10));
t44 = sin(qJ(2));
t40 = cos(pkin(6));
t47 = cos(qJ(2));
t70 = t40 * t47;
t23 = t37 * t44 - t39 * t70;
t43 = sin(qJ(4));
t46 = cos(qJ(4));
t73 = t38 * t46;
t10 = -t23 * t43 + t39 * t73;
t71 = t40 * t44;
t24 = t37 * t47 + t39 * t71;
t45 = cos(qJ(5));
t21 = t24 * t45;
t42 = sin(qJ(5));
t55 = t37 * t71 - t39 * t47;
t76 = t55 * t45;
t25 = t37 * t70 + t39 * t44;
t8 = t25 * t43 + t37 * t73;
t64 = t43 * t47;
t28 = -t38 * t64 + t40 * t46;
t63 = t44 * t45;
t86 = g(3) * (-t28 * t42 + t38 * t63);
t87 = -g(1) * (-t8 * t42 - t76) - g(2) * (t10 * t42 + t21) - t86;
t74 = t38 * t44;
t85 = g(3) * (-t28 * t45 - t42 * t74);
t66 = t43 * t44;
t84 = (-t42 * t66 + t45 * t47) * t88;
t83 = (t42 * t47 + t43 * t63) * t88;
t60 = t46 * t47;
t27 = t38 * t60 + t40 * t43;
t82 = g(3) * t27;
t81 = t37 * pkin(2);
t80 = t39 * pkin(2);
t19 = t23 * t42;
t78 = t24 * t42;
t77 = t55 * t42;
t75 = t38 * t43;
t72 = t38 * t47;
t48 = pkin(2) + pkin(8);
t69 = t40 * t48;
t68 = t42 * t25;
t67 = t42 * t43;
t65 = t43 * t45;
t62 = t44 * t46;
t61 = t45 * t25;
t35 = t37 * qJ(3);
t36 = t39 * qJ(3);
t59 = t40 * t36;
t54 = t40 * t60 - t75;
t58 = g(1) * (t54 * t37 + t39 * t62) - g(2) * (-t37 * t62 + t54 * t39);
t7 = -t25 * t46 + t37 * t75;
t9 = t23 * t46 + t39 * t75;
t57 = -g(1) * t7 + g(2) * t9;
t34 = t45 * pkin(5) + pkin(4);
t41 = -qJ(6) - pkin(9);
t56 = t43 * t34 + t46 * t41;
t53 = t40 * t64 + t73;
t51 = -t57 + t82;
t50 = g(1) * t8 - g(2) * t10 + g(3) * t28;
t2 = -g(1) * t25 - g(2) * t23 + g(3) * t72;
t49 = -g(1) * t55 + g(2) * t24 + g(3) * t74;
t33 = qJ(3) * t74;
t32 = t40 * t35;
t20 = t23 * t45;
t17 = t42 * t82;
t14 = t24 * t43;
t13 = t55 * t43;
t5 = -t37 * t66 + t53 * t39;
t3 = t53 * t37 + t39 * t66;
t1 = t49 * t46;
t4 = [-g(3), 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, -t2, t49, t2, -t49, -g(1) * (-(t32 + t80) * t44 + (-t40 * t81 + t36) * t47) - g(2) * (-(-t59 + t81) * t44 + (t40 * t80 + t35) * t47) - g(3) * (pkin(2) * t72 + t33), 0, 0, 0, 0, 0, -t49 * t43, -t1, 0, 0, 0, 0, 0, -g(1) * (-t13 * t45 - t68) - g(2) * (t14 * t45 - t19) - t83, -g(1) * (t13 * t42 - t61) - g(2) * (-t14 * t42 - t20) - t84, -g(1) * (-t55 * t65 - t68) - g(2) * (t24 * t65 - t19) - t83, -g(1) * (t55 * t67 - t61) - g(2) * (-t24 * t67 - t20) - t84, t1, -g(1) * (-pkin(5) * t68 - (t39 * t48 + t32) * t44 + (-t37 * t69 + t36) * t47 - t56 * t55) - g(2) * (-pkin(5) * t19 - (t37 * t48 - t59) * t44 + (t39 * t69 + t35) * t47 + t56 * t24) + (-t33 - ((pkin(5) * t42 + t48) * t47 + t56 * t44) * t38) * g(3); 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, t50, 0, 0, 0, 0, 0, (-t58 + t82) * t45, t58 * t42 - t17, t51 * t45, t57 * t42 - t17, -t50, -g(1) * (-t7 * t34 - t8 * t41) - g(2) * (t10 * t41 + t9 * t34) - g(3) * (-t27 * t34 - t28 * t41); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t3 * t42 - t76) - g(2) * (t5 * t42 + t21) - t86, -g(1) * (-t3 * t45 + t77) - g(2) * (t5 * t45 - t78) - t85, t87, -g(1) * (-t8 * t45 + t77) - g(2) * (t10 * t45 - t78) - t85, 0, t87 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51;];
taug_reg = t4;
