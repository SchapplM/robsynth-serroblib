% Calculate minimal parameter regressor of gravitation load for
% S6RPRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% taug_reg [6x29]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRR7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t13 = qJ(3) + pkin(10) + qJ(5);
t11 = sin(t13);
t12 = cos(t13);
t21 = sin(qJ(1));
t24 = cos(qJ(1));
t9 = g(1) * t21 - g(2) * t24;
t34 = -g(3) * t11 + t9 * t12;
t20 = sin(qJ(3));
t33 = pkin(3) * t20;
t31 = g(3) * t12;
t19 = sin(qJ(6));
t30 = t21 * t19;
t22 = cos(qJ(6));
t29 = t21 * t22;
t28 = t24 * t19;
t27 = t24 * t22;
t26 = t24 * pkin(1) + t21 * qJ(2);
t10 = g(1) * t24 + g(2) * t21;
t23 = cos(qJ(3));
t25 = g(3) * t20 - t9 * t23;
t18 = -qJ(4) - pkin(7);
t15 = t24 * qJ(2);
t8 = t11 * t27 - t30;
t7 = t11 * t28 + t29;
t6 = t11 * t29 + t28;
t5 = -t11 * t30 + t27;
t3 = t9 * t11 + t31;
t2 = t34 * t22;
t1 = t34 * t19;
t4 = [0, t9, t10, -t9, -t10, -g(1) * (-t21 * pkin(1) + t15) - g(2) * t26, 0, 0, 0, 0, 0, -t10 * t20, -t10 * t23, t9, -g(1) * (t24 * t33 + t15 + (-pkin(1) + t18) * t21) - g(2) * (-t24 * t18 + t21 * t33 + t26) 0, 0, 0, 0, 0, -t10 * t11, -t10 * t12, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t6, g(1) * t7 - g(2) * t5; 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, g(3) * t23 + t9 * t20, 0, t25 * pkin(3), 0, 0, 0, 0, 0, -t34, t3, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, t3, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 - g(2) * t7 + t19 * t31, g(1) * t6 - g(2) * t8 + t22 * t31;];
taug_reg  = t4;
