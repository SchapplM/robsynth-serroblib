% Calculate minimal parameter regressor of gravitation load for
% S5RPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% taug_reg [5x23]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRR8_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR8_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t12 = qJ(4) + qJ(5);
t10 = sin(t12);
t17 = sin(qJ(3));
t18 = sin(qJ(1));
t19 = cos(qJ(3));
t20 = cos(qJ(1));
t3 = -t18 * t17 - t20 * t19;
t4 = t20 * t17 - t18 * t19;
t16 = g(1) * t4 - g(2) * t3;
t24 = t16 * t10;
t11 = cos(t12);
t23 = t16 * t11;
t13 = sin(qJ(4));
t22 = t16 * t13;
t14 = cos(qJ(4));
t21 = t16 * t14;
t15 = g(1) * t3 + g(2) * t4;
t6 = g(1) * t20 + g(2) * t18;
t5 = g(1) * t18 - g(2) * t20;
t2 = -g(3) * t10 - t15 * t11;
t1 = g(3) * t11 - t15 * t10;
t7 = [0, t5, t6, t5, -t6, -g(1) * (-t18 * pkin(1) + t20 * qJ(2)) - g(2) * (t20 * pkin(1) + t18 * qJ(2)), 0, -t16, t15, 0, 0, 0, 0, 0, -t21, t22, 0, 0, 0, 0, 0, -t23, t24; 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t16, -t15, 0, 0, 0, 0, 0, t21, -t22, 0, 0, 0, 0, 0, t23, -t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t14 - t15 * t13, -g(3) * t13 - t15 * t14, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t7;
