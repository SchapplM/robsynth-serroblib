% Calculate minimal parameter regressor of gravitation load for
% S5RRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [5x31]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRR12_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR12_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR12_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR12_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t49 = sin(qJ(2));
t52 = cos(qJ(2));
t53 = cos(qJ(1));
t66 = cos(pkin(5));
t59 = t53 * t66;
t77 = sin(qJ(1));
t33 = t77 * t49 - t52 * t59;
t34 = t49 * t59 + t77 * t52;
t48 = sin(qJ(3));
t45 = cos(pkin(6));
t78 = cos(qJ(3));
t62 = t45 * t78;
t43 = sin(pkin(6));
t44 = sin(pkin(5));
t72 = t44 * t53;
t64 = t43 * t72;
t12 = t33 * t62 + t34 * t48 + t78 * t64;
t71 = t45 * t48;
t13 = -t33 * t71 + t34 * t78 - t48 * t64;
t25 = -t33 * t43 + t45 * t72;
t47 = sin(qJ(4));
t51 = cos(qJ(4));
t4 = t13 * t51 - t25 * t47;
t46 = sin(qJ(5));
t50 = cos(qJ(5));
t83 = -t12 * t50 + t4 * t46;
t82 = t12 * t46 + t4 * t50;
t79 = t13 * t47 + t25 * t51;
t76 = t43 * t47;
t75 = t43 * t51;
t74 = t44 * t49;
t73 = t44 * t52;
t70 = t46 * t51;
t69 = t48 * t49;
t68 = t48 * t52;
t67 = t50 * t51;
t65 = t43 * t74;
t63 = t44 * t77;
t61 = t78 * t49;
t60 = t78 * t52;
t58 = t66 * t43;
t57 = t43 * t63;
t56 = t66 * t77;
t24 = t48 * t58 + (t45 * t68 + t61) * t44;
t32 = -t43 * t73 + t66 * t45;
t35 = -t53 * t49 - t52 * t56;
t36 = -t49 * t56 + t53 * t52;
t17 = t36 * t78 + (t35 * t45 + t57) * t48;
t27 = -t35 * t43 + t45 * t63;
t6 = -t17 * t47 + t27 * t51;
t55 = g(1) * t6 - g(2) * t79 + g(3) * (-t24 * t47 + t32 * t51);
t16 = -t35 * t62 + t36 * t48 - t78 * t57;
t23 = -t78 * t58 + (-t45 * t60 + t69) * t44;
t54 = g(1) * t16 + g(2) * t12 + g(3) * t23;
t31 = (-t45 * t69 + t60) * t44;
t30 = (t45 * t61 + t68) * t44;
t22 = t31 * t51 + t47 * t65;
t21 = t35 * t78 - t36 * t71;
t20 = t35 * t48 + t36 * t62;
t19 = -t33 * t78 - t34 * t71;
t18 = -t33 * t48 + t34 * t62;
t11 = t24 * t51 + t32 * t47;
t9 = t21 * t51 + t36 * t76;
t8 = t19 * t51 + t34 * t76;
t7 = t17 * t51 + t27 * t47;
t2 = t16 * t46 + t7 * t50;
t1 = t16 * t50 - t7 * t46;
t3 = [0, g(1) * t77 - g(2) * t53, g(1) * t53 + g(2) * t77, 0, 0, 0, 0, 0, g(1) * t34 - g(2) * t36, -g(1) * t33 - g(2) * t35, 0, 0, 0, 0, 0, g(1) * t13 - g(2) * t17, -g(1) * t12 + g(2) * t16, 0, 0, 0, 0, 0, g(1) * t4 - g(2) * t7, -g(1) * t79 - g(2) * t6, 0, 0, 0, 0, 0, g(1) * t82 - g(2) * t2, -g(1) * t83 - g(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t35 + g(2) * t33 - g(3) * t73, g(1) * t36 + g(2) * t34 + g(3) * t74, 0, 0, 0, 0, 0, -g(1) * t21 - g(2) * t19 - g(3) * t31, g(1) * t20 + g(2) * t18 + g(3) * t30, 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t8 - g(3) * t22, -g(1) * (-t21 * t47 + t36 * t75) - g(2) * (-t19 * t47 + t34 * t75) - g(3) * (-t31 * t47 + t51 * t65), 0, 0, 0, 0, 0, -g(1) * (t20 * t46 + t9 * t50) - g(2) * (t18 * t46 + t8 * t50) - g(3) * (t22 * t50 + t30 * t46), -g(1) * (t20 * t50 - t9 * t46) - g(2) * (t18 * t50 - t8 * t46) - g(3) * (-t22 * t46 + t30 * t50); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, g(1) * t17 + g(2) * t13 + g(3) * t24, 0, 0, 0, 0, 0, t54 * t51, -t54 * t47, 0, 0, 0, 0, 0, -g(1) * (-t16 * t67 + t17 * t46) - g(2) * (-t12 * t67 + t13 * t46) - g(3) * (-t23 * t67 + t24 * t46), -g(1) * (t16 * t70 + t17 * t50) - g(2) * (t12 * t70 + t13 * t50) - g(3) * (t23 * t70 + t24 * t50); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, g(1) * t7 + g(2) * t4 + g(3) * t11, 0, 0, 0, 0, 0, -t55 * t50, t55 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 + g(2) * t83 - g(3) * (-t11 * t46 + t23 * t50), g(1) * t2 + g(2) * t82 - g(3) * (-t11 * t50 - t23 * t46);];
taug_reg = t3;
