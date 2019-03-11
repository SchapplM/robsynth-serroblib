% Calculate minimal parameter regressor of gravitation load for
% S6RPRRPP7
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
% Datum: 2019-03-09 04:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRPP7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t70 = -pkin(1) - pkin(7);
t69 = -pkin(4) - pkin(5);
t39 = sin(qJ(1));
t68 = g(1) * t39;
t41 = cos(qJ(3));
t67 = g(2) * t41;
t42 = cos(qJ(1));
t66 = g(2) * t42;
t38 = sin(qJ(3));
t65 = g(3) * t38;
t33 = t41 * pkin(8);
t37 = sin(qJ(4));
t64 = t37 * t41;
t63 = t38 * t39;
t62 = t39 * t37;
t40 = cos(qJ(4));
t61 = t39 * t40;
t60 = t39 * t41;
t59 = t40 * t41;
t58 = t42 * t37;
t57 = t42 * t40;
t56 = -pkin(8) + qJ(6);
t55 = t42 * pkin(1) + t39 * qJ(2);
t54 = t41 * qJ(6);
t53 = t37 * t60;
t52 = t39 * t59;
t51 = -qJ(5) * t37 - pkin(3);
t14 = t38 * t62 - t57;
t15 = t38 * t61 + t58;
t50 = -t14 * pkin(4) + t15 * qJ(5);
t16 = t38 * t58 + t61;
t17 = t38 * t57 - t62;
t49 = t16 * pkin(4) - t17 * qJ(5);
t48 = pkin(3) * t60 + pkin(4) * t52 + pkin(8) * t63 + qJ(5) * t53;
t47 = g(1) * t16 + g(2) * t14;
t23 = g(1) * t42 + g(2) * t39;
t22 = -t66 + t68;
t46 = -pkin(4) * t40 + t51;
t45 = pkin(3) * t63 + t15 * pkin(4) + t42 * pkin(7) + t14 * qJ(5) + t55;
t2 = g(1) * t14 - g(2) * t16 + g(3) * t64;
t44 = g(1) * t15 - g(2) * t17 + g(3) * t59;
t32 = t42 * qJ(2);
t43 = t17 * pkin(4) + t16 * qJ(5) + t32 + (t38 * pkin(3) - t33) * t42;
t9 = -t22 * t41 + t65;
t24 = qJ(5) * t59;
t18 = t23 * t41;
t8 = g(3) * t41 + t22 * t38;
t7 = t9 * t40;
t6 = g(1) * t53 - t37 * t65 - t58 * t67;
t5 = -g(1) * t17 - g(2) * t15;
t1 = [0, t22, t23, -t22, -t23, -g(1) * (-t39 * pkin(1) + t32) - g(2) * t55, 0, 0, 0, 0, 0, -t23 * t38, -t18, 0, 0, 0, 0, 0, t5, t47, t5, t18, -t47, -g(1) * (t70 * t39 + t43) - g(2) * (-pkin(8) * t60 + t45) t5, -t47, -t18, -g(1) * (t17 * pkin(5) + t42 * t54 + t43) - g(2) * (t15 * pkin(5) + t45) + (-g(1) * t70 - t56 * t67) * t39; 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t8, 0, 0, 0, 0, 0, t7, t6, t7, -t8, -t6, -g(1) * t48 - g(3) * t33 - t46 * t65 - (-pkin(8) * t38 + t46 * t41) * t66, t7, -t6, t8, -g(1) * (pkin(5) * t52 + t48) - g(3) * (t33 - t54) + (qJ(6) * t68 - g(3) * (-pkin(5) * t40 + t46)) * t38 - (t56 * t38 + (t69 * t40 + t51) * t41) * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t44, t2, 0, -t44, -g(1) * t50 - g(2) * t49 - g(3) * (-pkin(4) * t64 + t24) t2, -t44, 0, -g(1) * (-t14 * pkin(5) + t50) - g(2) * (t16 * pkin(5) + t49) - g(3) * (t69 * t64 + t24); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9;];
taug_reg  = t1;
