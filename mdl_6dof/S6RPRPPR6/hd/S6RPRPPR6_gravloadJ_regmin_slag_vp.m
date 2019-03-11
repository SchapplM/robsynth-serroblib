% Calculate minimal parameter regressor of gravitation load for
% S6RPRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
% 
% Output:
% taug_reg [6x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPPR6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t24 = sin(qJ(1));
t26 = cos(qJ(1));
t49 = -g(1) * t24 + g(2) * t26;
t19 = qJ(3) + pkin(9);
t11 = sin(t19);
t13 = cos(t19);
t29 = -g(3) * t11 - t13 * t49;
t45 = g(3) * t13;
t23 = sin(qJ(3));
t44 = t23 * pkin(3);
t18 = pkin(10) + qJ(6);
t10 = sin(t18);
t43 = t24 * t10;
t12 = cos(t18);
t42 = t24 * t12;
t20 = sin(pkin(10));
t41 = t24 * t20;
t21 = cos(pkin(10));
t40 = t24 * t21;
t39 = t26 * t10;
t38 = t26 * t12;
t37 = t26 * t20;
t36 = t26 * t21;
t35 = t26 * pkin(1) + t24 * qJ(2);
t34 = -t24 * pkin(1) + t26 * qJ(2);
t6 = g(1) * t26 + g(2) * t24;
t22 = -qJ(4) - pkin(7);
t33 = t24 * t22 + t26 * t44 + t34;
t32 = t11 * pkin(4) - t13 * qJ(5);
t31 = -t26 * t22 + t24 * t44 + t35;
t25 = cos(qJ(3));
t27 = g(3) * t23 + t25 * t49;
t4 = t11 * t38 - t43;
t3 = t11 * t39 + t42;
t2 = t11 * t42 + t39;
t1 = -t11 * t43 + t38;
t5 = [0, -t49, t6, t49, -t6, -g(1) * t34 - g(2) * t35, 0, 0, 0, 0, 0, -t6 * t23, -t6 * t25, -t49, -g(1) * t33 - g(2) * t31, -g(1) * (t11 * t36 - t41) - g(2) * (t11 * t40 + t37) -g(1) * (-t11 * t37 - t40) - g(2) * (-t11 * t41 + t36) t6 * t13, -g(1) * (t32 * t26 + t33) - g(2) * (t32 * t24 + t31) 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t2, g(1) * t3 - g(2) * t1; 0, 0, 0, 0, 0, t49, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0, 0, 0, t49, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, g(3) * t25 - t23 * t49, 0, t27 * pkin(3), -t29 * t21, t29 * t20, t11 * t49 - t45, -g(3) * (-t32 - t44) + t49 * (pkin(3) * t25 + pkin(4) * t13 + qJ(5) * t11) 0, 0, 0, 0, 0, -t29 * t12, t29 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, 0, 0, 0, -t6, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 - g(2) * t3 + t10 * t45, g(1) * t2 - g(2) * t4 + t12 * t45;];
taug_reg  = t5;
