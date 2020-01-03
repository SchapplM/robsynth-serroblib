% Calculate minimal parameter regressor of gravitation load for
% S5RPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% taug_reg [5x19]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRR8_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR8_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t16 = sin(qJ(5));
t18 = cos(qJ(1));
t24 = pkin(8) + qJ(4);
t22 = sin(t24);
t23 = cos(t24);
t26 = sin(qJ(1));
t1 = -t18 * t23 - t26 * t22;
t2 = t18 * t22 - t26 * t23;
t21 = g(1) * t2 - g(2) * t1;
t28 = t21 * t16;
t17 = cos(qJ(5));
t27 = t21 * t17;
t25 = t18 * pkin(1) + t26 * qJ(2);
t20 = g(1) * t1 + g(2) * t2;
t19 = -t26 * pkin(1) + t18 * qJ(2);
t15 = cos(pkin(8));
t14 = sin(pkin(8));
t6 = g(1) * t18 + g(2) * t26;
t5 = g(1) * t26 - g(2) * t18;
t4 = t26 * t14 + t18 * t15;
t3 = t18 * t14 - t26 * t15;
t7 = [0, t5, t6, t5, -t6, -g(1) * t19 - g(2) * t25, -g(1) * t3 - g(2) * t4, -g(1) * t4 + g(2) * t3, -g(1) * (-t26 * pkin(2) + t19) - g(2) * (t18 * pkin(2) + t25), 0, -t21, t20, 0, 0, 0, 0, 0, -t27, t28; 0, 0, 0, 0, 0, -t5, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t20, 0, 0, 0, 0, 0, t27, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t17 - t20 * t16, -g(3) * t16 - t20 * t17;];
taug_reg = t7;
