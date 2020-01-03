% Calculate minimal parameter regressor of gravitation load for
% S5RRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% taug_reg [5x23]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPPR7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t18 = sin(qJ(1));
t21 = cos(qJ(1));
t7 = g(1) * t21 + g(2) * t18;
t14 = qJ(2) + pkin(8);
t11 = cos(t14);
t30 = g(3) * t11;
t16 = sin(qJ(5));
t29 = t18 * t16;
t19 = cos(qJ(5));
t28 = t18 * t19;
t27 = t21 * t16;
t26 = t21 * t19;
t6 = g(1) * t18 - g(2) * t21;
t10 = sin(t14);
t25 = t11 * pkin(3) + t10 * qJ(4);
t23 = -g(3) * t10 - t7 * t11;
t17 = sin(qJ(2));
t20 = cos(qJ(2));
t22 = -g(3) * t20 + t7 * t17;
t15 = -qJ(3) - pkin(6);
t12 = t20 * pkin(2);
t9 = t12 + pkin(1);
t8 = t21 * t9;
t5 = -t10 * t29 + t26;
t4 = t10 * t28 + t27;
t3 = t10 * t27 + t28;
t2 = t10 * t26 - t29;
t1 = -t7 * t10 + t30;
t13 = [0, t6, t7, 0, 0, 0, 0, 0, t6 * t20, -t6 * t17, -t7, -g(1) * (-t21 * t15 - t18 * t9) - g(2) * (-t18 * t15 + t8), -t7, -t6 * t11, t6 * t10, -g(2) * t8 + (g(1) * t15 - g(2) * t25) * t21 + (-g(1) * (-t25 - t9) + g(2) * t15) * t18, 0, 0, 0, 0, 0, -g(1) * t5 - g(2) * t3, g(1) * t4 - g(2) * t2; 0, 0, 0, 0, 0, 0, 0, 0, t22, g(3) * t17 + t7 * t20, 0, t22 * pkin(2), 0, t1, t23, -g(3) * (t12 + t25) + t7 * (pkin(2) * t17 + pkin(3) * t10 - qJ(4) * t11), 0, 0, 0, 0, 0, t23 * t16, t23 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, 0, 0, 0, -t6, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t4 + t19 * t30, g(1) * t3 - g(2) * t5 - t16 * t30;];
taug_reg = t13;
