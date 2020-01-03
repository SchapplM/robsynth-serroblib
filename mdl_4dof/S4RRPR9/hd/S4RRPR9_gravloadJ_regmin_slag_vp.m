% Calculate minimal parameter regressor of gravitation load for
% S4RRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% taug_reg [4x21]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRPR9_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR9_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR9_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR9_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t14 = sin(qJ(1));
t16 = cos(qJ(1));
t22 = g(1) * t16 + g(2) * t14;
t13 = sin(qJ(2));
t15 = cos(qJ(2));
t5 = -g(3) * t15 + t22 * t13;
t28 = g(3) * t13;
t26 = t14 * t15;
t25 = t15 * t16;
t11 = sin(pkin(7));
t24 = t16 * t11;
t12 = cos(pkin(7));
t23 = t16 * t12;
t21 = g(1) * t14 - g(2) * t16;
t20 = t15 * pkin(2) + t13 * qJ(3);
t18 = pkin(1) + t20;
t10 = pkin(7) + qJ(4);
t9 = cos(t10);
t8 = sin(t10);
t7 = t21 * t13;
t6 = t22 * t15 + t28;
t4 = t14 * t8 + t9 * t25;
t3 = t14 * t9 - t8 * t25;
t2 = t16 * t8 - t9 * t26;
t1 = t16 * t9 + t8 * t26;
t17 = [0, t21, t22, 0, 0, 0, 0, 0, t21 * t15, -t7, -g(1) * (-t12 * t26 + t24) - g(2) * (t14 * t11 + t15 * t23), -g(1) * (t11 * t26 + t23) - g(2) * (t14 * t12 - t15 * t24), t7, (-g(1) * pkin(5) - g(2) * t18) * t16 + (-g(2) * pkin(5) + g(1) * t18) * t14, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t4, -g(1) * t1 - g(2) * t3; 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, t5 * t12, -t5 * t11, -t6, -g(3) * t20 + t22 * (pkin(2) * t13 - qJ(3) * t15), 0, 0, 0, 0, 0, t5 * t9, -t5 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 + g(2) * t1 + t8 * t28, g(1) * t4 - g(2) * t2 + t9 * t28;];
taug_reg = t17;
