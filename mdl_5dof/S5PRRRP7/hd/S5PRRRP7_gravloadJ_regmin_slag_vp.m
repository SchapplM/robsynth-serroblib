% Calculate minimal parameter regressor of gravitation load for
% S5PRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% taug_reg [5x22]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:46
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRRP7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:45:23
% EndTime: 2021-01-15 16:45:26
% DurationCPUTime: 0.40s
% Computational Cost: add. (225->93), mult. (591->153), div. (0->0), fcn. (723->10), ass. (0->65)
t38 = sin(pkin(9));
t40 = cos(pkin(9));
t48 = cos(qJ(2));
t41 = cos(pkin(5));
t45 = sin(qJ(2));
t64 = t41 * t45;
t26 = t38 * t48 + t40 * t64;
t44 = sin(qJ(3));
t39 = sin(pkin(5));
t47 = cos(qJ(3));
t56 = t47 * t39;
t3 = t26 * t44 + t40 * t56;
t24 = t38 * t64 - t40 * t48;
t5 = t24 * t44 + t38 * t56;
t82 = g(1) * t5 - g(2) * t3;
t72 = g(3) * t39;
t59 = t44 * t45;
t63 = t41 * t47;
t29 = t39 * t59 - t63;
t73 = g(3) * t29;
t81 = -t73 + t82;
t62 = t41 * t48;
t25 = t38 * t45 - t40 * t62;
t46 = cos(qJ(4));
t18 = t25 * t46;
t60 = t44 * t39;
t2 = t24 * t47 - t38 * t60;
t27 = t38 * t62 + t40 * t45;
t21 = t27 * t46;
t4 = t26 * t47 - t40 * t60;
t43 = sin(qJ(4));
t30 = t41 * t44 + t45 * t56;
t66 = t39 * t48;
t77 = g(3) * (-t30 * t43 - t46 * t66);
t80 = -g(1) * (t2 * t43 + t21) - g(2) * (-t4 * t43 + t18) - t77;
t76 = g(3) * (-t30 * t46 + t43 * t66);
t55 = t47 * t48;
t75 = (-t43 * t55 + t45 * t46) * t72;
t74 = (t43 * t45 + t46 * t55) * t72;
t22 = t24 * t43;
t69 = t25 * t43;
t19 = t26 * t43;
t68 = t27 * t43;
t67 = t38 * t41;
t65 = t40 * t41;
t61 = t43 * t47;
t58 = t44 * t48;
t57 = t46 * t47;
t28 = -t41 * t59 - t56;
t54 = -g(1) * (t38 * t28 + t40 * t58) - g(2) * (-t40 * t28 + t38 * t58);
t37 = t46 * pkin(4) + pkin(3);
t42 = -qJ(5) - pkin(8);
t52 = -t47 * t37 + t44 * t42;
t51 = g(1) * t2 - g(2) * t4 - g(3) * t30;
t49 = g(1) * t27 + g(2) * t25 - g(3) * t66;
t31 = t45 * t63 - t60;
t23 = t24 * t46;
t20 = t26 * t46;
t17 = t43 * t73;
t14 = t27 * t47;
t13 = t25 * t47;
t10 = -t38 * t31 + t40 * t55;
t8 = t40 * t31 + t38 * t55;
t1 = t49 * t44;
t6 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, t49, -g(1) * t24 + g(2) * t26 + t45 * t72, 0, 0, 0, 0, 0, t49 * t47, -t1, 0, 0, 0, 0, 0, -g(1) * (-t14 * t46 - t22) - g(2) * (-t13 * t46 + t19) - t74, -g(1) * (t14 * t43 - t23) - g(2) * (t13 * t43 + t20) - t75, -g(1) * (-t27 * t57 - t22) - g(2) * (-t25 * t57 + t19) - t74, -g(1) * (t27 * t61 - t23) - g(2) * (t25 * t61 + t20) - t75, t1, -g(1) * (-pkin(4) * t22 - (t40 * pkin(2) + pkin(7) * t67) * t45 + (-pkin(2) * t67 + t40 * pkin(7)) * t48 + t52 * t27) - g(2) * (pkin(4) * t19 - (t38 * pkin(2) - pkin(7) * t65) * t45 + (pkin(2) * t65 + t38 * pkin(7)) * t48 + t52 * t25) - ((pkin(4) * t43 + pkin(7)) * t45 + (pkin(2) - t52) * t48) * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, -t51, 0, 0, 0, 0, 0, (-t54 + t73) * t46, t54 * t43 - t17, -t81 * t46, t82 * t43 - t17, t51, -g(1) * (t2 * t42 + t5 * t37) - g(2) * (-t3 * t37 - t4 * t42) - g(3) * (-t29 * t37 - t30 * t42); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t10 * t43 + t21) - g(2) * (-t8 * t43 + t18) - t77, -g(1) * (-t10 * t46 - t68) - g(2) * (-t8 * t46 - t69) - t76, t80, -g(1) * (t2 * t46 - t68) - g(2) * (-t4 * t46 - t69) - t76, 0, t80 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81;];
taug_reg = t6;
