% Calculate minimal parameter regressor of gravitation load for
% S5PPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% taug_reg [5x14]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PPRRP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t8 = sin(pkin(7));
t9 = cos(pkin(7));
t15 = g(1) * t9 + g(2) * t8;
t7 = pkin(8) + qJ(3);
t5 = sin(t7);
t6 = cos(t7);
t1 = -g(3) * t6 + t15 * t5;
t21 = g(3) * t5;
t11 = sin(qJ(4));
t19 = t8 * t11;
t12 = cos(qJ(4));
t18 = t8 * t12;
t17 = t9 * t11;
t16 = t9 * t12;
t13 = -g(1) * (-t6 * t17 + t18) - g(2) * (-t6 * t19 - t16) + t11 * t21;
t10 = -qJ(5) - pkin(6);
t4 = t12 * pkin(4) + pkin(3);
t3 = -g(1) * t8 + g(2) * t9;
t2 = t15 * t6 + t21;
t14 = [-g(3), -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3; 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, t1 * t12, -t1 * t11, -t2, -g(3) * (-t5 * t10 + t6 * t4) + t15 * (t10 * t6 + t4 * t5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -g(1) * (-t6 * t16 - t19) - g(2) * (-t6 * t18 + t17) + t12 * t21, 0, t13 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t14;
