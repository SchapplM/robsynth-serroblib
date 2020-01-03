% Calculate minimal parameter regressor of gravitation load for
% S5RPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% taug_reg [5x22]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR14_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR14_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t17 = sin(qJ(1));
t20 = cos(qJ(1));
t5 = g(1) * t17 - g(2) * t20;
t13 = qJ(3) + pkin(8);
t7 = sin(t13);
t8 = cos(t13);
t30 = -g(3) * t7 + t5 * t8;
t28 = g(3) * t8;
t16 = sin(qJ(3));
t27 = pkin(3) * t16;
t26 = t20 * pkin(1) + t17 * qJ(2);
t15 = sin(qJ(5));
t25 = t17 * t15;
t18 = cos(qJ(5));
t24 = t17 * t18;
t23 = t20 * t15;
t22 = t20 * t18;
t6 = g(1) * t20 + g(2) * t17;
t19 = cos(qJ(3));
t21 = g(3) * t16 - t5 * t19;
t14 = -qJ(4) - pkin(6);
t10 = t20 * qJ(2);
t4 = t7 * t22 - t25;
t3 = t7 * t23 + t24;
t2 = t7 * t24 + t23;
t1 = -t7 * t25 + t22;
t9 = [0, t5, t6, -t5, -t6, -g(1) * (-t17 * pkin(1) + t10) - g(2) * t26, 0, 0, 0, 0, 0, -t6 * t16, -t6 * t19, t5, -g(1) * (t20 * t27 + t10 + (-pkin(1) + t14) * t17) - g(2) * (-t20 * t14 + t17 * t27 + t26), 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t2, g(1) * t3 - g(2) * t1; 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, g(3) * t19 + t5 * t16, 0, t21 * pkin(3), 0, 0, 0, 0, 0, -t30 * t18, t30 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 - g(2) * t3 + t15 * t28, g(1) * t2 - g(2) * t4 + t18 * t28;];
taug_reg = t9;
