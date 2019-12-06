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
% Datum: 2019-12-05 18:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
t20 = pkin(9) + qJ(4);
t21 = qJ(1) + qJ(2);
t18 = sin(t21);
t19 = cos(t21);
t27 = -t18 * pkin(2) + t19 * qJ(3);
t10 = g(2) * t19 + g(3) * t18;
t9 = g(2) * t18 - g(3) * t19;
t26 = -t19 * pkin(2) - t18 * qJ(3);
t25 = cos(qJ(1));
t24 = sin(qJ(1));
t17 = qJ(5) + t20;
t16 = cos(t20);
t15 = sin(t20);
t13 = cos(t17);
t12 = sin(t17);
t8 = t10 * cos(pkin(9));
t7 = t10 * sin(pkin(9));
t6 = t10 * t16;
t5 = t10 * t15;
t4 = t10 * t13;
t3 = t10 * t12;
t2 = -g(1) * t13 - t9 * t12;
t1 = g(1) * t12 - t9 * t13;
t11 = [0, g(2) * t25 + g(3) * t24, -g(2) * t24 + g(3) * t25, 0, t10, -t9, t8, -t7, t9, -g(2) * (-t25 * pkin(1) + t26) - g(3) * (-t24 * pkin(1) + t27), 0, 0, 0, 0, 0, t6, -t5, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, t10, -t9, t8, -t7, t9, -g(2) * t26 - g(3) * t27, 0, 0, 0, 0, 0, t6, -t5, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t16 - t9 * t15, g(1) * t15 - t9 * t16, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1;];
taug_reg = t11;
