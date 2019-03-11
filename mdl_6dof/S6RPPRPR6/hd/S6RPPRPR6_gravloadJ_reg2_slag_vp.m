% Calculate inertial parameters regressor of gravitation load for
% S6RPPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRPR6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t26 = sin(qJ(1));
t29 = cos(qJ(1));
t10 = g(1) * t29 + g(2) * t26;
t25 = sin(qJ(4));
t28 = cos(qJ(4));
t1 = g(3) * t28 + t10 * t25;
t50 = pkin(5) + pkin(7);
t49 = g(3) * t25;
t47 = t25 * pkin(4);
t46 = t25 * pkin(8);
t45 = t29 * pkin(7);
t44 = t26 * t28;
t43 = t28 * t29;
t24 = sin(qJ(6));
t42 = t29 * t24;
t27 = cos(qJ(6));
t41 = t29 * t27;
t40 = -pkin(1) - qJ(3);
t18 = t28 * qJ(5);
t20 = t29 * qJ(2);
t39 = t26 * t18 + t20;
t35 = qJ(5) * t25;
t38 = pkin(4) * t44 + t26 * t35;
t37 = pkin(4) * t43 + t29 * t35;
t36 = t29 * pkin(1) + t26 * qJ(2);
t34 = t29 * qJ(3) + t36;
t33 = t29 * t47 + t34;
t32 = t40 - t47;
t9 = g(1) * t26 - g(2) * t29;
t31 = t40 * t26 + t20;
t8 = t9 * t28;
t7 = t9 * t25;
t6 = t24 * t44 - t41;
t5 = t27 * t44 + t42;
t4 = t26 * t27 + t28 * t42;
t3 = t26 * t24 - t28 * t41;
t2 = t10 * t28 - t49;
t11 = [0, 0, 0, 0, 0, 0, t9, t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, -t10, -g(1) * (-t26 * pkin(1) + t20) - g(2) * t36, 0, 0, 0, 0, 0, 0, 0, -t10, t9, -g(1) * t31 - g(2) * t34, 0, 0, 0, 0, 0, 0, t7, t8, t10, -g(1) * (t31 - t45) - g(2) * (-t26 * pkin(7) + t34) 0, 0, 0, 0, 0, 0, t10, -t7, -t8, -g(1) * (t39 - t45) - g(2) * (-t29 * t18 + t33) + (g(2) * pkin(7) - g(1) * t32) * t26, 0, 0, 0, 0, 0, 0, -g(1) * t6 + g(2) * t4, -g(1) * t5 - g(2) * t3, t7, -g(1) * t39 - g(2) * t33 + (g(1) * t50 - g(2) * (-t18 + t46)) * t29 + (-g(1) * (t32 - t46) + g(2) * t50) * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1, -g(1) * t37 - g(2) * t38 - g(3) * (t18 - t47) 0, 0, 0, 0, 0, 0, -t1 * t24, -t1 * t27, -t2, -g(1) * (pkin(8) * t43 + t37) - g(2) * (pkin(8) * t44 + t38) - g(3) * (t18 + (-pkin(4) - pkin(8)) * t25); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 + g(2) * t5 - t27 * t49, -g(1) * t4 - g(2) * t6 + t24 * t49, 0, 0;];
taug_reg  = t11;
