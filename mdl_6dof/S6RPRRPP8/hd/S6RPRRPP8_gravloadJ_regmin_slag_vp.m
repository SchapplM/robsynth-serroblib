% Calculate minimal parameter regressor of gravitation load for
% S6RPRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
% 
% Output:
% taug_reg [6x28]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRPP8_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP8_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t44 = cos(qJ(3));
t45 = cos(qJ(1));
t63 = t44 * t45;
t41 = sin(qJ(3));
t71 = g(3) * t41;
t76 = -g(2) * t63 - t71;
t75 = -pkin(1) - pkin(7);
t74 = -pkin(5) - pkin(8);
t42 = sin(qJ(1));
t73 = g(1) * t42;
t72 = g(2) * t45;
t40 = sin(qJ(4));
t70 = t40 * t44;
t69 = t41 * t42;
t68 = t41 * t45;
t67 = t42 * t40;
t43 = cos(qJ(4));
t66 = t42 * t43;
t65 = t42 * t44;
t64 = t43 * t44;
t62 = t45 * t40;
t61 = t45 * t43;
t60 = -pkin(4) - qJ(6);
t59 = t45 * pkin(1) + t42 * qJ(2);
t57 = t40 * t65;
t56 = t42 * t64;
t55 = -qJ(5) * t40 - pkin(3);
t13 = t41 * t67 - t61;
t14 = t41 * t66 + t62;
t54 = -t13 * pkin(4) + t14 * qJ(5);
t15 = t41 * t62 + t66;
t16 = t41 * t61 - t67;
t53 = t15 * pkin(4) - t16 * qJ(5);
t52 = pkin(3) * t65 + pkin(4) * t56 + pkin(8) * t69 + qJ(5) * t57;
t51 = g(1) * t15 + g(2) * t13;
t50 = g(1) * t16 + g(2) * t14;
t23 = g(1) * t45 + g(2) * t42;
t22 = -t72 + t73;
t49 = -pkin(4) * t43 + t55;
t48 = pkin(3) * t69 + t14 * pkin(4) + t45 * pkin(7) + t13 * qJ(5) + t59;
t2 = g(1) * t13 - g(2) * t15 + g(3) * t70;
t47 = g(1) * t14 - g(2) * t16 + g(3) * t64;
t35 = t45 * qJ(2);
t46 = pkin(3) * t68 + t16 * pkin(4) - pkin(8) * t63 + t15 * qJ(5) + t35;
t36 = t44 * pkin(8);
t24 = qJ(5) * t64;
t17 = t23 * t44;
t8 = g(1) * t69 - g(2) * t68 + g(3) * t44;
t7 = g(1) * t56 + t76 * t43;
t6 = g(1) * t57 + t76 * t40;
t1 = [0, t22, t23, -t22, -t23, -g(1) * (-t42 * pkin(1) + t35) - g(2) * t59, 0, 0, 0, 0, 0, -t23 * t41, -t17, 0, 0, 0, 0, 0, -t50, t51, t17, t50, -t51, -g(1) * (t75 * t42 + t46) - g(2) * (-pkin(8) * t65 + t48) t17, -t51, -t50, -g(1) * (-pkin(5) * t63 + t16 * qJ(6) + t46) - g(2) * (t14 * qJ(6) + t48) + (-g(2) * t74 * t44 - g(1) * t75) * t42; 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22 * t44 + t71, t8, 0, 0, 0, 0, 0, -t7, t6, -t8, t7, -t6, -g(1) * t52 - g(3) * t36 - t49 * t71 - (-pkin(8) * t41 + t49 * t44) * t72, -t8, -t6, -t7, -g(1) * (qJ(6) * t56 + t52) - g(3) * (t44 * pkin(5) + t36) + (-pkin(5) * t73 - g(3) * (-qJ(6) * t43 + t49)) * t41 - (t74 * t41 + (t60 * t43 + t55) * t44) * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t47, 0, -t2, -t47, -g(1) * t54 - g(2) * t53 - g(3) * (-pkin(4) * t70 + t24) 0, -t47, t2, -g(1) * (-t13 * qJ(6) + t54) - g(2) * (t15 * qJ(6) + t53) - g(3) * (t60 * t70 + t24); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47;];
taug_reg  = t1;
