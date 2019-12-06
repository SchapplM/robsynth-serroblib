% Calculate minimal parameter regressor of gravitation load for
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% taug_reg [5x17]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t15 = qJ(1) + qJ(2);
t14 = qJ(3) + t15;
t9 = pkin(9) + t14;
t7 = sin(t9);
t8 = cos(t9);
t23 = g(2) * t8 + g(3) * t7;
t22 = -g(2) * t7 + g(3) * t8;
t10 = sin(t14);
t12 = sin(t15);
t21 = -pkin(2) * t12 - pkin(3) * t10;
t11 = cos(t14);
t13 = cos(t15);
t20 = -pkin(2) * t13 - pkin(3) * t11;
t4 = g(2) * t11 + g(3) * t10;
t19 = cos(qJ(1));
t18 = cos(qJ(5));
t17 = sin(qJ(1));
t16 = sin(qJ(5));
t6 = g(2) * t13 + g(3) * t12;
t5 = -g(2) * t12 + g(3) * t13;
t3 = -g(2) * t10 + g(3) * t11;
t2 = t23 * t18;
t1 = t23 * t16;
t24 = [0, g(2) * t19 + g(3) * t17, -g(2) * t17 + g(3) * t19, 0, t6, t5, 0, t4, t3, -g(2) * (-t19 * pkin(1) + t20) - g(3) * (-t17 * pkin(1) + t21), 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, t6, t5, 0, t4, t3, -g(2) * t20 - g(3) * t21, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, t4, t3, t4 * pkin(3), 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t18 + t22 * t16, g(1) * t16 + t22 * t18;];
taug_reg = t24;
