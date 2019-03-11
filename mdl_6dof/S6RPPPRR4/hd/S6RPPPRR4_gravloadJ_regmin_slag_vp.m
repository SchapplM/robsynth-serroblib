% Calculate minimal parameter regressor of gravitation load for
% S6RPPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
% 
% Output:
% taug_reg [6x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPPRR4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t17 = sin(qJ(5));
t19 = cos(qJ(5));
t28 = sin(pkin(9));
t29 = cos(pkin(9));
t33 = sin(qJ(1));
t34 = cos(qJ(1));
t4 = -t33 * t28 - t34 * t29;
t5 = t34 * t28 - t33 * t29;
t26 = g(1) * t5 - g(2) * t4;
t37 = g(3) * t17 + t26 * t19;
t35 = g(3) * t19;
t16 = sin(qJ(6));
t32 = t16 * t17;
t18 = cos(qJ(6));
t31 = t17 * t18;
t30 = t34 * pkin(1) + t33 * qJ(2);
t27 = t34 * pkin(2) + t30;
t25 = g(1) * t4 + g(2) * t5;
t24 = -t33 * pkin(1) + t34 * qJ(2);
t23 = t5 * t16 + t4 * t31;
t22 = -t5 * t18 + t4 * t32;
t21 = -t33 * pkin(2) + t24;
t7 = g(1) * t34 + g(2) * t33;
t6 = g(1) * t33 - g(2) * t34;
t2 = -t4 * t16 + t5 * t31;
t1 = -t4 * t18 - t5 * t32;
t3 = [0, t6, t7, t6, -t7, -g(1) * t24 - g(2) * t30, -t26, t25, -g(1) * t21 - g(2) * t27, t26, -t25, -g(1) * (t5 * pkin(3) + t4 * qJ(4) + t21) - g(2) * (-t4 * pkin(3) + t5 * qJ(4) + t27) 0, 0, 0, 0, 0, -t25 * t17, -t25 * t19, 0, 0, 0, 0, 0, -g(1) * t23 - g(2) * t2, g(1) * t22 - g(2) * t1; 0, 0, 0, 0, 0, -t6, 0, 0, -t6, 0, 0, -t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, t26 * t17 - t35, 0, 0, 0, 0, 0, -t37 * t18, t37 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 - g(2) * t22 - t16 * t35, g(1) * t2 - g(2) * t23 - t18 * t35;];
taug_reg  = t3;
