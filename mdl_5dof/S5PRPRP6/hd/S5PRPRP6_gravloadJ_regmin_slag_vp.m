% Calculate minimal parameter regressor of gravitation load for
% S5PRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
% 
% Output:
% taug_reg [5x18]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRPRP6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t16 = sin(pkin(7));
t17 = cos(pkin(7));
t35 = -g(1) * t17 - g(2) * t16;
t19 = sin(qJ(2));
t21 = cos(qJ(2));
t5 = g(3) * t19 - t35 * t21;
t34 = pkin(2) * t19;
t30 = g(3) * t21;
t18 = sin(qJ(4));
t29 = t18 * t19;
t20 = cos(qJ(4));
t28 = t19 * t20;
t27 = t21 * pkin(2) + t19 * qJ(3);
t26 = qJ(3) * t21;
t24 = pkin(4) * t18 - qJ(5) * t20;
t7 = t16 * t20 + t17 * t29;
t9 = -t16 * t29 + t17 * t20;
t23 = -g(1) * t7 + g(2) * t9 + t18 * t30;
t6 = t16 * t18 - t17 * t28;
t8 = t16 * t28 + t17 * t18;
t1 = g(1) * t6 - g(2) * t8 + t20 * t30;
t12 = t17 * t26;
t11 = t16 * t26;
t4 = -t35 * t19 - t30;
t3 = t5 * t20;
t2 = t5 * t18;
t10 = [-g(3), 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, t4, t5, -t4, -t5, -g(1) * (-t17 * t34 + t12) - g(2) * (-t16 * t34 + t11) - g(3) * t27, 0, 0, 0, 0, 0, -t2, -t3, -t2, t4, t3, -g(1) * t12 - g(2) * t11 - g(3) * (t21 * pkin(6) + t24 * t19 + t27) + t35 * (t24 * t21 + (-pkin(2) - pkin(6)) * t19); 0, 0, 0, 0, 0, 0, -t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t23, t1, 0, t23, -g(1) * (-t6 * pkin(4) + t7 * qJ(5)) - g(2) * (t8 * pkin(4) - t9 * qJ(5)) - (-pkin(4) * t20 - qJ(5) * t18) * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t10;
