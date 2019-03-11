% Calculate minimal parameter regressor of gravitation load for
% S6RPPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
% 
% Output:
% taug_reg [6x27]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRPR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t20 = sin(qJ(4));
t22 = cos(qJ(4));
t25 = -t20 * pkin(4) + t22 * qJ(5);
t21 = sin(qJ(1));
t23 = cos(qJ(1));
t9 = g(1) * t23 + g(2) * t21;
t6 = -g(3) * t20 + t9 * t22;
t40 = g(3) * t22;
t17 = pkin(9) + qJ(6);
t10 = sin(t17);
t38 = t21 * t10;
t11 = cos(t17);
t37 = t21 * t11;
t18 = sin(pkin(9));
t36 = t21 * t18;
t19 = cos(pkin(9));
t35 = t21 * t19;
t34 = t23 * t10;
t33 = t23 * t11;
t32 = t23 * t18;
t31 = t23 * t19;
t30 = -pkin(1) - qJ(3);
t29 = t23 * pkin(1) + t21 * qJ(2);
t27 = t23 * qJ(3) + t29;
t8 = g(1) * t21 - g(2) * t23;
t14 = t23 * qJ(2);
t7 = t8 * t22;
t5 = t9 * t20 + t40;
t4 = t20 * t33 - t38;
t3 = -t20 * t34 - t37;
t2 = -t20 * t37 - t34;
t1 = t20 * t38 - t33;
t12 = [0, t8, t9, -t8, -t9, -g(1) * (-t21 * pkin(1) + t14) - g(2) * t29, -t9, t8, -g(1) * (t30 * t21 + t14) - g(2) * t27, 0, 0, 0, 0, 0, t8 * t20, t7, -g(1) * (-t20 * t35 - t32) - g(2) * (t20 * t31 - t36) -g(1) * (t20 * t36 - t31) - g(2) * (-t20 * t32 - t35) -t7, -g(1) * (-t23 * pkin(7) + t14) - g(2) * (-t25 * t23 + t27) + (-g(1) * (t25 + t30) + g(2) * pkin(7)) * t21, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t4, -g(1) * t1 - g(2) * t3; 0, 0, 0, 0, 0, -t8, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5, -t6 * t19, t6 * t18, -t5, -g(3) * t25 - t9 * (pkin(4) * t22 + qJ(5) * t20) 0, 0, 0, 0, 0, -t6 * t11, t6 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 + g(2) * t1 + t10 * t40, g(1) * t4 - g(2) * t2 + t11 * t40;];
taug_reg  = t12;
