% Calculate minimal parameter regressor of gravitation load for
% S5RRPPR3
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
% taug_reg [5x17]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPPR3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR3_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR3_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t17 = qJ(1) + qJ(2);
t14 = sin(t17);
t26 = pkin(2) * t14;
t19 = sin(qJ(1));
t25 = t19 * pkin(1);
t13 = pkin(8) + t17;
t10 = sin(t13);
t11 = cos(t13);
t15 = cos(t17);
t12 = pkin(2) * t15;
t24 = t11 * pkin(3) + t10 * qJ(4) + t12;
t4 = -g(1) * t11 - g(2) * t10;
t23 = g(1) * t10 - g(2) * t11;
t5 = g(1) * t14 - g(2) * t15;
t22 = -t10 * pkin(3) + t11 * qJ(4) - t26;
t21 = cos(qJ(1));
t20 = cos(qJ(5));
t18 = sin(qJ(5));
t16 = t21 * pkin(1);
t6 = g(1) * t15 + g(2) * t14;
t2 = t4 * t20;
t1 = t4 * t18;
t3 = [0, g(1) * t19 - g(2) * t21, g(1) * t21 + g(2) * t19, 0, t5, t6, -g(1) * (-t25 - t26) - g(2) * (t12 + t16), -t23, t4, -g(1) * (t22 - t25) - g(2) * (t16 + t24), 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, t5, t6, t5 * pkin(2), -t23, t4, -g(1) * t22 - g(2) * t24, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, -g(3), 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t18 - t23 * t20, g(3) * t20 + t23 * t18;];
taug_reg = t3;
