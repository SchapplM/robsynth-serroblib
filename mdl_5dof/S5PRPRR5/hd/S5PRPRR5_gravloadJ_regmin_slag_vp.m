% Calculate minimal parameter regressor of gravitation load for
% S5PRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% taug_reg [5x22]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:55
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRPRR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR5_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR5_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t12 = sin(pkin(8));
t14 = cos(pkin(8));
t19 = g(1) * t14 + g(2) * t12;
t15 = sin(qJ(2));
t16 = cos(qJ(2));
t3 = -g(3) * t16 + t19 * t15;
t23 = g(3) * t15;
t21 = t12 * t16;
t20 = t14 * t16;
t10 = pkin(9) + qJ(4);
t9 = qJ(5) + t10;
t8 = cos(t10);
t7 = sin(t10);
t6 = cos(t9);
t5 = sin(t9);
t4 = t19 * t16 + t23;
t2 = -g(1) * (-t12 * t5 - t6 * t20) - g(2) * (t14 * t5 - t6 * t21) + t6 * t23;
t1 = -g(1) * (t12 * t6 - t5 * t20) - g(2) * (-t14 * t6 - t5 * t21) + t5 * t23;
t11 = [-g(3), 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t3, t4, t3 * cos(pkin(9)), -t3 * sin(pkin(9)), -t4, -g(3) * (t16 * pkin(2) + t15 * qJ(3)) + t19 * (pkin(2) * t15 - qJ(3) * t16), 0, 0, 0, 0, 0, t3 * t8, -t3 * t7, 0, 0, 0, 0, 0, t3 * t6, -t3 * t5; 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t12 * t8 - t7 * t20) - g(2) * (-t14 * t8 - t7 * t21) + t7 * t23, -g(1) * (-t12 * t7 - t8 * t20) - g(2) * (t14 * t7 - t8 * t21) + t8 * t23, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t11;
