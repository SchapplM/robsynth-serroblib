% Calculate minimal parameter regressor of gravitation load for
% S5PRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% 
% Output:
% taug_reg [5x21]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRPRR7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t10 = cos(pkin(8));
t9 = sin(pkin(8));
t17 = g(1) * t10 + g(2) * t9;
t12 = sin(qJ(2));
t14 = cos(qJ(2));
t4 = g(3) * t12 + t17 * t14;
t22 = g(3) * t14;
t21 = t12 * t9;
t20 = t10 * t12;
t11 = sin(qJ(4));
t19 = t11 * t12;
t13 = cos(qJ(4));
t18 = t12 * t13;
t8 = qJ(4) + qJ(5);
t6 = cos(t8);
t5 = sin(t8);
t3 = t17 * t12 - t22;
t2 = -g(1) * (-t5 * t20 - t9 * t6) - g(2) * (t10 * t6 - t5 * t21) - t5 * t22;
t1 = -g(1) * (t6 * t20 - t9 * t5) - g(2) * (t10 * t5 + t6 * t21) + t6 * t22;
t7 = [-g(3), 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t3, t4, -t3, -t4, -g(3) * (t14 * pkin(2) + t12 * qJ(3)) + t17 * (pkin(2) * t12 - qJ(3) * t14), 0, 0, 0, 0, 0, -t4 * t11, -t4 * t13, 0, 0, 0, 0, 0, -t4 * t5, -t4 * t6; 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t10 * t18 - t9 * t11) - g(2) * (t10 * t11 + t9 * t18) + t13 * t22, -g(1) * (-t10 * t19 - t9 * t13) - g(2) * (t10 * t13 - t9 * t19) - t11 * t22, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t7;
