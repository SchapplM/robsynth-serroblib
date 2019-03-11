% Calculate minimal parameter regressor of gravitation load for
% S6RPPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
% 
% Output:
% taug_reg [6x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPPRR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t15 = sin(qJ(5));
t18 = cos(qJ(5));
t13 = cos(pkin(9));
t16 = sin(qJ(1));
t19 = cos(qJ(1));
t27 = sin(pkin(9));
t3 = t19 * t13 - t16 * t27;
t4 = t16 * t13 + t19 * t27;
t25 = g(1) * t4 - g(2) * t3;
t20 = -g(3) * t18 + t25 * t15;
t32 = g(3) * t15;
t30 = t19 * pkin(1) + t16 * qJ(2);
t14 = sin(qJ(6));
t29 = t14 * t18;
t17 = cos(qJ(6));
t28 = t17 * t18;
t26 = t19 * qJ(3) + t30;
t24 = g(1) * t3 + g(2) * t4;
t10 = t19 * qJ(2);
t23 = t10 + (-pkin(1) - qJ(3)) * t16;
t22 = t4 * t14 + t3 * t28;
t21 = -t4 * t17 + t3 * t29;
t6 = g(1) * t19 + g(2) * t16;
t5 = g(1) * t16 - g(2) * t19;
t2 = -t3 * t14 + t4 * t28;
t1 = -t3 * t17 - t4 * t29;
t7 = [0, t5, t6, -t5, -t6, -g(1) * (-t16 * pkin(1) + t10) - g(2) * t30, -t6, t5, -g(1) * t23 - g(2) * t26, -t24, t25, -g(1) * (t19 * pkin(3) + t23) - g(2) * (t16 * pkin(3) + t26) 0, 0, 0, 0, 0, -t24 * t18, t24 * t15, 0, 0, 0, 0, 0, -g(1) * t22 - g(2) * t2, g(1) * t21 - g(2) * t1; 0, 0, 0, 0, 0, -t5, 0, 0, -t5, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t6, 0, 0, -t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, t25 * t18 + t32, 0, 0, 0, 0, 0, t20 * t17, -t20 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 - g(2) * t21 + t14 * t32, g(1) * t2 - g(2) * t22 + t17 * t32;];
taug_reg  = t7;
