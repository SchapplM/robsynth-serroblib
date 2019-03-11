% Calculate minimal parameter regressor of gravitation load for
% S6RPRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
% 
% Output:
% taug_reg [6x24]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRPP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t36 = cos(qJ(4));
t23 = t36 * pkin(4) + pkin(3);
t32 = -qJ(5) - pkin(8);
t34 = sin(qJ(3));
t37 = cos(qJ(3));
t48 = t37 * t23 - t34 * t32;
t64 = -pkin(2) - t48;
t31 = qJ(1) + pkin(9);
t25 = sin(t31);
t27 = cos(t31);
t47 = g(1) * t27 + g(2) * t25;
t5 = -g(3) * t37 + t47 * t34;
t61 = g(3) * t34;
t33 = sin(qJ(4));
t59 = t25 * t33;
t58 = t25 * t36;
t57 = t25 * t37;
t56 = t27 * t33;
t55 = t27 * t36;
t54 = t27 * t37;
t53 = t32 * t37;
t52 = t33 * t37;
t50 = t36 * t37;
t49 = t27 * t52;
t46 = g(1) * t25 - g(2) * t27;
t35 = sin(qJ(1));
t38 = cos(qJ(1));
t45 = g(1) * t35 - g(2) * t38;
t30 = qJ(4) + pkin(10);
t24 = sin(t30);
t26 = cos(t30);
t44 = pkin(5) * t26 + qJ(6) * t24;
t7 = t25 * t52 + t55;
t1 = t24 * t57 + t27 * t26;
t3 = t24 * t54 - t25 * t26;
t42 = g(1) * t3 + g(2) * t1 + t24 * t61;
t41 = t38 * pkin(1) + pkin(4) * t59 + t25 * pkin(7) - t64 * t27;
t6 = t47 * t37 + t61;
t39 = -t35 * pkin(1) + pkin(4) * t56 + t27 * pkin(7) + t64 * t25;
t17 = pkin(4) * t58;
t11 = t46 * t34;
t10 = t27 * t50 + t59;
t9 = -t49 + t58;
t8 = -t25 * t50 + t56;
t4 = t25 * t24 + t26 * t54;
t2 = -t27 * t24 + t26 * t57;
t12 = [0, t45, g(1) * t38 + g(2) * t35, t45 * pkin(1), 0, 0, 0, 0, 0, t46 * t37, -t11, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t10, -g(1) * t7 - g(2) * t9, t11, -g(1) * t39 - g(2) * t41, g(1) * t2 - g(2) * t4, t11, g(1) * t1 - g(2) * t3, -g(1) * (-t2 * pkin(5) - t1 * qJ(6) + t39) - g(2) * (t4 * pkin(5) + t3 * qJ(6) + t41); 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t5 * t36, -t5 * t33, -t6, -g(3) * t48 + t47 * (t23 * t34 + t53) t5 * t26, -t6, t5 * t24, -g(3) * (t44 * t37 + t48) + t47 * (t53 - (-t23 - t44) * t34); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 + g(2) * t7 + t33 * t61, g(1) * t10 - g(2) * t8 + t36 * t61, 0, -g(1) * t17 + (g(2) * t55 + t33 * t6) * pkin(4), t42, 0, -g(1) * t4 - g(2) * t2 - t26 * t61, -g(1) * (-pkin(4) * t49 - t3 * pkin(5) + t4 * qJ(6) + t17) - g(2) * (-t7 * pkin(4) - t1 * pkin(5) + t2 * qJ(6)) - (-pkin(4) * t33 - pkin(5) * t24 + qJ(6) * t26) * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42;];
taug_reg  = t12;
