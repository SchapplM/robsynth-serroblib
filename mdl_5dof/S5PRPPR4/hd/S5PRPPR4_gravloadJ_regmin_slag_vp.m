% Calculate minimal parameter regressor of gravitation load for
% S5PRPPR4
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
% taug_reg [5x19]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRPPR4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR4_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR4_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t15 = pkin(7) + qJ(2);
t13 = sin(t15);
t24 = g(1) * t13;
t14 = cos(t15);
t23 = t14 * pkin(2) + t13 * qJ(3);
t8 = g(1) * t14 + g(2) * t13;
t7 = -g(2) * t14 + t24;
t16 = sin(pkin(8));
t17 = cos(pkin(8));
t22 = pkin(3) * t17 + qJ(4) * t16;
t18 = sin(qJ(5));
t19 = cos(qJ(5));
t21 = t16 * t19 - t17 * t18;
t20 = t16 * t18 + t17 * t19;
t10 = t14 * qJ(3);
t6 = t7 * t17;
t5 = t7 * t16;
t4 = t20 * t14;
t3 = t21 * t14;
t2 = t20 * t13;
t1 = t21 * t13;
t9 = [-g(3), 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, t7, t8, t6, -t5, -t8, -g(1) * (-t13 * pkin(2) + t10) - g(2) * t23, t6, -t8, t5, -g(1) * t10 - g(2) * (t22 * t14 + t23) - (-pkin(2) - t22) * t24, 0, 0, 0, 0, 0, g(1) * t2 - g(2) * t4, g(1) * t1 - g(2) * t3; 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t17 - t8 * t16, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 - g(2) * t1 + g(3) * t20, g(1) * t4 + g(2) * t2 + g(3) * t21;];
taug_reg = t9;
