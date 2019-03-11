% Calculate minimal parameter regressor of gravitation load for
% S6RPRRPP4
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
% taug_reg [6x27]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRPP4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t39 = cos(qJ(4));
t25 = t39 * pkin(4) + pkin(3);
t31 = pkin(9) + qJ(3);
t28 = cos(t31);
t13 = t28 * t25;
t26 = sin(t31);
t35 = -qJ(5) - pkin(8);
t62 = -t26 * t35 + t13;
t34 = cos(pkin(9));
t61 = -t34 * pkin(2) - pkin(1) - t62;
t38 = sin(qJ(1));
t40 = cos(qJ(1));
t17 = g(1) * t40 + g(2) * t38;
t1 = -g(3) * t28 + t17 * t26;
t58 = g(3) * t26;
t32 = qJ(4) + pkin(10);
t27 = sin(t32);
t55 = t38 * t27;
t29 = cos(t32);
t54 = t38 * t29;
t37 = sin(qJ(4));
t53 = t38 * t37;
t52 = t38 * t39;
t51 = t40 * t27;
t50 = t40 * t29;
t49 = t40 * t37;
t48 = t40 * t39;
t47 = t28 * t49;
t16 = g(1) * t38 - g(2) * t40;
t46 = pkin(5) * t29 + qJ(6) * t27;
t8 = t28 * t53 + t48;
t3 = t28 * t55 + t50;
t5 = t28 * t51 - t54;
t43 = g(1) * t5 + g(2) * t3 + t27 * t58;
t36 = -pkin(7) - qJ(2);
t42 = pkin(4) * t53 - t38 * t36 - t61 * t40;
t2 = t17 * t28 + t58;
t41 = pkin(4) * t49 - t40 * t36 + t61 * t38;
t21 = pkin(4) * t52;
t11 = t28 * t48 + t53;
t10 = -t47 + t52;
t9 = -t28 * t52 + t49;
t7 = t16 * t26;
t6 = t28 * t50 + t55;
t4 = t28 * t54 - t51;
t12 = [0, t16, t17, t16 * t34, -t16 * sin(pkin(9)) -t17, -g(1) * (-t38 * pkin(1) + t40 * qJ(2)) - g(2) * (t40 * pkin(1) + t38 * qJ(2)) 0, 0, 0, 0, 0, t16 * t28, -t7, 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t11, -g(1) * t8 - g(2) * t10, t7, -g(1) * t41 - g(2) * t42, g(1) * t4 - g(2) * t6, t7, g(1) * t3 - g(2) * t5, -g(1) * (-t4 * pkin(5) - t3 * qJ(6) + t41) - g(2) * (t6 * pkin(5) + t5 * qJ(6) + t42); 0, 0, 0, 0, 0, 0, -t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, 0, 0, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, t1 * t39, -t1 * t37, -t2, -g(3) * t62 + t17 * (t25 * t26 + t28 * t35) t1 * t29, -t2, t1 * t27, -g(3) * t13 + (-g(3) * t46 + t17 * t35) * t28 + (g(3) * t35 + t17 * (t25 + t46)) * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t10 + g(2) * t8 + t37 * t58, g(1) * t11 - g(2) * t9 + t39 * t58, 0, -g(1) * t21 + (g(2) * t48 + t2 * t37) * pkin(4), t43, 0, -g(1) * t6 - g(2) * t4 - t29 * t58, -g(1) * (-pkin(4) * t47 - t5 * pkin(5) + t6 * qJ(6) + t21) - g(2) * (-t8 * pkin(4) - t3 * pkin(5) + t4 * qJ(6)) - (-pkin(4) * t37 - pkin(5) * t27 + qJ(6) * t29) * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43;];
taug_reg  = t12;
