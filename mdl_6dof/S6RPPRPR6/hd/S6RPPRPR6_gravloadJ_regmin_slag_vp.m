% Calculate minimal parameter regressor of gravitation load for
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
% taug_reg [6x27]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRPR6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t17 = sin(qJ(4));
t20 = cos(qJ(4));
t23 = -t17 * pkin(4) + t20 * qJ(5);
t18 = sin(qJ(1));
t21 = cos(qJ(1));
t10 = g(1) * t21 + g(2) * t18;
t1 = g(3) * t20 + t10 * t17;
t34 = g(3) * t17;
t31 = t18 * t20;
t16 = sin(qJ(6));
t30 = t21 * t16;
t19 = cos(qJ(6));
t29 = t21 * t19;
t28 = -pkin(1) - qJ(3);
t27 = t21 * pkin(1) + t18 * qJ(2);
t25 = t21 * qJ(3) + t27;
t9 = g(1) * t18 - g(2) * t21;
t13 = t21 * qJ(2);
t8 = t9 * t20;
t7 = t9 * t17;
t6 = t16 * t31 - t29;
t5 = t19 * t31 + t30;
t4 = t18 * t19 + t20 * t30;
t3 = t18 * t16 - t20 * t29;
t2 = t10 * t20 - t34;
t11 = [0, t9, t10, -t9, -t10, -g(1) * (-t18 * pkin(1) + t13) - g(2) * t27, -t10, t9, -g(1) * (t28 * t18 + t13) - g(2) * t25, 0, 0, 0, 0, 0, t7, t8, t10, -t7, -t8, -g(1) * (-t21 * pkin(7) + t13) - g(2) * (-t23 * t21 + t25) + (-g(1) * (t23 + t28) + g(2) * pkin(7)) * t18, 0, 0, 0, 0, 0, -g(1) * t6 + g(2) * t4, -g(1) * t5 - g(2) * t3; 0, 0, 0, 0, 0, -t9, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, 0, t2, -t1, -g(3) * t23 - t10 * (pkin(4) * t20 + qJ(5) * t17) 0, 0, 0, 0, 0, -t1 * t16, -t1 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 + g(2) * t5 - t19 * t34, -g(1) * t4 - g(2) * t6 + t16 * t34;];
taug_reg  = t11;
