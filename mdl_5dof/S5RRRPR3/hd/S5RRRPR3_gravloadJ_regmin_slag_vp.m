% Calculate minimal parameter regressor of gravitation load for
% S5RRRPR3
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
% taug_reg [5x22]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPR3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t21 = cos(qJ(3));
t13 = t21 * pkin(3) + pkin(2);
t17 = qJ(1) + qJ(2);
t15 = sin(t17);
t16 = cos(t17);
t18 = -qJ(4) - pkin(7);
t25 = -t16 * t13 + t15 * t18;
t8 = g(2) * t16 + g(3) * t15;
t7 = g(2) * t15 - g(3) * t16;
t24 = -t15 * t13 - t16 * t18;
t19 = sin(qJ(3));
t23 = -g(1) * t21 - t7 * t19;
t22 = cos(qJ(1));
t20 = sin(qJ(1));
t14 = qJ(3) + pkin(9) + qJ(5);
t11 = cos(t14);
t10 = sin(t14);
t6 = t8 * t21;
t5 = t8 * t19;
t4 = t8 * t11;
t3 = t8 * t10;
t2 = -g(1) * t11 - t7 * t10;
t1 = g(1) * t10 - t7 * t11;
t9 = [0, g(2) * t22 + g(3) * t20, -g(2) * t20 + g(3) * t22, 0, t8, -t7, 0, 0, 0, 0, 0, t6, -t5, t7, -g(2) * (-t22 * pkin(1) + t25) - g(3) * (-t20 * pkin(1) + t24), 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, t8, -t7, 0, 0, 0, 0, 0, t6, -t5, t7, -g(2) * t25 - g(3) * t24, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, g(1) * t19 - t7 * t21, 0, t23 * pkin(3), 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1;];
taug_reg = t9;
