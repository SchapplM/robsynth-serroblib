% Calculate minimal parameter regressor of gravitation load for
% S5RRPRR4
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
% taug_reg [5x21]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t17 = qJ(1) + qJ(2);
t11 = pkin(9) + t17;
t10 = cos(t11);
t9 = sin(t11);
t23 = -g(2) * t9 + g(3) * t10;
t22 = g(2) * t10 + g(3) * t9;
t13 = sin(t17);
t15 = cos(t17);
t8 = g(2) * t15 + g(3) * t13;
t21 = cos(qJ(1));
t20 = cos(qJ(4));
t19 = sin(qJ(1));
t18 = sin(qJ(4));
t16 = qJ(4) + qJ(5);
t14 = cos(t16);
t12 = sin(t16);
t7 = -g(2) * t13 + g(3) * t15;
t6 = t22 * t20;
t5 = t22 * t18;
t4 = t22 * t14;
t3 = t22 * t12;
t2 = -g(1) * t14 + t23 * t12;
t1 = g(1) * t12 + t23 * t14;
t24 = [0, g(2) * t21 + g(3) * t19, -g(2) * t19 + g(3) * t21, 0, t8, t7, -g(2) * (-t21 * pkin(1) - pkin(2) * t15) - g(3) * (-t19 * pkin(1) - pkin(2) * t13), 0, 0, 0, 0, 0, t6, -t5, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, t8, t7, t8 * pkin(2), 0, 0, 0, 0, 0, t6, -t5, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t20 + t23 * t18, g(1) * t18 + t23 * t20, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1;];
taug_reg = t24;
