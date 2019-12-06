% Calculate minimal parameter regressor of gravitation load for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% taug_reg [5x16]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t13 = qJ(1) + qJ(2);
t11 = sin(t13);
t26 = pkin(2) * t11;
t12 = cos(t13);
t25 = pkin(2) * t12;
t16 = sin(qJ(1));
t24 = t16 * pkin(1);
t18 = cos(qJ(1));
t23 = t18 * pkin(1);
t10 = pkin(8) + t13;
t7 = sin(t10);
t8 = cos(t10);
t22 = g(2) * t8 + g(3) * t7;
t3 = g(2) * t7 - g(3) * t8;
t5 = g(2) * t12 + g(3) * t11;
t14 = -qJ(5) - pkin(7);
t17 = cos(qJ(4));
t9 = t17 * pkin(4) + pkin(3);
t21 = t7 * t14 - t8 * t9 - t25;
t20 = -t8 * t14 - t7 * t9 - t26;
t15 = sin(qJ(4));
t19 = -g(1) * t17 - t3 * t15;
t4 = -g(2) * t11 + g(3) * t12;
t2 = t22 * t17;
t1 = t22 * t15;
t6 = [0, g(2) * t18 + g(3) * t16, -g(2) * t16 + g(3) * t18, 0, t5, t4, -g(2) * (-t23 - t25) - g(3) * (-t24 - t26), 0, 0, 0, 0, 0, t2, -t1, t3, -g(2) * (t21 - t23) - g(3) * (t20 - t24); 0, 0, 0, 0, t5, t4, t5 * pkin(2), 0, 0, 0, 0, 0, t2, -t1, t3, -g(2) * t21 - g(3) * t20; 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, -g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, g(1) * t15 - t3 * t17, 0, t19 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22;];
taug_reg = t6;
