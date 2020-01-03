% Calculate minimal parameter regressor of gravitation load for
% S5RRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% taug_reg [5x24]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t23 = qJ(1) + qJ(2);
t20 = sin(t23);
t21 = cos(t23);
t29 = t21 * pkin(2) + t20 * qJ(3);
t22 = pkin(9) + qJ(4);
t28 = t20 * pkin(2) - t21 * qJ(3);
t10 = g(2) * t21 + g(3) * t20;
t9 = g(2) * t20 - g(3) * t21;
t27 = cos(qJ(1));
t26 = sin(qJ(1));
t19 = qJ(5) + t22;
t18 = cos(t22);
t17 = sin(t22);
t13 = cos(t19);
t12 = sin(t19);
t8 = t10 * cos(pkin(9));
t7 = t10 * sin(pkin(9));
t6 = t10 * t18;
t5 = t10 * t17;
t4 = t10 * t13;
t3 = t10 * t12;
t2 = -g(1) * t13 + t9 * t12;
t1 = g(1) * t12 + t9 * t13;
t11 = [0, -g(2) * t27 - g(3) * t26, g(2) * t26 - g(3) * t27, 0, -t10, t9, -t8, t7, -t9, -g(2) * (t27 * pkin(1) + t29) - g(3) * (t26 * pkin(1) + t28), 0, 0, 0, 0, 0, -t6, t5, 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, -t10, t9, -t8, t7, -t9, -g(2) * t29 - g(3) * t28, 0, 0, 0, 0, 0, -t6, t5, 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t18 + t9 * t17, g(1) * t17 + t9 * t18, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1;];
taug_reg = t11;
