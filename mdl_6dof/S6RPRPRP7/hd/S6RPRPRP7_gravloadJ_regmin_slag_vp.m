% Calculate minimal parameter regressor of gravitation load for
% S6RPRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% 
% Output:
% taug_reg [6x24]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRP7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t19 = -qJ(4) - pkin(7);
t22 = sin(qJ(1));
t25 = cos(qJ(1));
t21 = sin(qJ(3));
t38 = t21 * pkin(3);
t48 = t22 * t19 + t25 * t38;
t47 = -g(1) * t22 + g(2) * t25;
t17 = qJ(3) + pkin(9);
t11 = sin(t17);
t23 = cos(qJ(5));
t34 = t25 * t23;
t20 = sin(qJ(5));
t37 = t22 * t20;
t1 = -t11 * t37 + t34;
t35 = t25 * t20;
t36 = t22 * t23;
t3 = t11 * t35 + t36;
t12 = cos(t17);
t39 = g(3) * t12;
t46 = -g(1) * t1 - g(2) * t3 + t20 * t39;
t27 = -g(3) * t11 - t12 * t47;
t43 = pkin(5) * t20;
t33 = t25 * pkin(1) + t22 * qJ(2);
t31 = t22 * t38 + t33;
t14 = t25 * qJ(2);
t30 = -t22 * pkin(1) + t14;
t6 = g(1) * t25 + g(2) * t22;
t10 = t23 * pkin(5) + pkin(4);
t18 = -qJ(6) - pkin(8);
t29 = t11 * t10 + t12 * t18;
t24 = cos(qJ(3));
t26 = g(3) * t21 + t24 * t47;
t4 = t11 * t34 - t37;
t2 = t11 * t36 + t35;
t5 = [0, -t47, t6, t47, -t6, -g(1) * t30 - g(2) * t33, 0, 0, 0, 0, 0, -t6 * t21, -t6 * t24, -t47, -g(1) * (t30 + t48) - g(2) * (-t25 * t19 + t31) 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t2, g(1) * t3 - g(2) * t1, t6 * t12, -g(1) * (t14 + t48) - g(2) * t31 + (-g(1) * t29 - g(2) * (-t19 + t43)) * t25 + (-g(1) * (-pkin(1) - t43) - g(2) * t29) * t22; 0, 0, 0, 0, 0, t47, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, 0, 0, 0, 0, 0, 0, 0, t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, g(3) * t24 - t21 * t47, 0, t26 * pkin(3), 0, 0, 0, 0, 0, -t27 * t23, t27 * t20, t11 * t47 - t39, -g(3) * (-t29 - t38) + t47 * (pkin(3) * t24 + t10 * t12 - t11 * t18); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, 0, 0, 0, 0, 0, 0, 0, 0, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, g(1) * t2 - g(2) * t4 + t23 * t39, 0, t46 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27;];
taug_reg  = t5;
