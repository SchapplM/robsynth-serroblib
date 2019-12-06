% Calculate minimal parameter regressor of gravitation load for
% S5PRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% taug_reg [5x18]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRPR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR1_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR1_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t18 = pkin(8) + qJ(2);
t16 = qJ(3) + t18;
t10 = sin(t16);
t11 = cos(t16);
t22 = t11 * pkin(3) + t10 * qJ(4);
t21 = -t10 * pkin(3) + t11 * qJ(4);
t6 = g(1) * t11 + g(2) * t10;
t5 = g(1) * t10 - g(2) * t11;
t17 = pkin(9) + qJ(5);
t15 = cos(t18);
t14 = cos(t17);
t13 = sin(t18);
t12 = sin(t17);
t4 = t5 * cos(pkin(9));
t3 = t5 * sin(pkin(9));
t2 = t5 * t14;
t1 = t5 * t12;
t7 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, g(1) * t13 - g(2) * t15, g(1) * t15 + g(2) * t13, 0, t5, t6, t4, -t3, -t6, -g(1) * (-pkin(2) * t13 + t21) - g(2) * (pkin(2) * t15 + t22), 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, t5, t6, t4, -t3, -t6, -g(1) * t21 - g(2) * t22, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t14 + t6 * t12, g(3) * t12 + t6 * t14;];
taug_reg = t7;
