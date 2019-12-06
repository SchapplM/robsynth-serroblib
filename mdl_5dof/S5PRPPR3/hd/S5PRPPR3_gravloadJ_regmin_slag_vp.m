% Calculate minimal parameter regressor of gravitation load for
% S5PRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
% 
% Output:
% taug_reg [5x15]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRPPR3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t6 = sin(pkin(7));
t7 = cos(pkin(7));
t15 = g(1) * t7 + g(2) * t6;
t5 = qJ(2) + pkin(8);
t4 = cos(t5);
t20 = g(3) * t4;
t8 = sin(qJ(5));
t19 = t6 * t8;
t18 = t7 * t8;
t10 = cos(qJ(5));
t17 = t6 * t10;
t16 = t7 * t10;
t3 = sin(t5);
t13 = -g(3) * t3 - t15 * t4;
t11 = cos(qJ(2));
t9 = sin(qJ(2));
t12 = -g(3) * t11 + t15 * t9;
t2 = -g(1) * t6 + g(2) * t7;
t1 = -t15 * t3 + t20;
t14 = [-g(3), 0, 0, 0, -g(3), 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, t12, g(3) * t9 + t15 * t11, t12 * pkin(2), t1, t13, -g(3) * (t11 * pkin(2) + t4 * pkin(3) + t3 * qJ(4)) + t15 * (pkin(2) * t9 + pkin(3) * t3 - qJ(4) * t4), 0, 0, 0, 0, 0, t13 * t8, t13 * t10; 0, 0, 0, 0, t2, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t3 * t16 - t19) - g(2) * (t3 * t17 + t18) + t10 * t20, -g(1) * (-t3 * t18 - t17) - g(2) * (-t3 * t19 + t16) - t8 * t20;];
taug_reg = t14;
