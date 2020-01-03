% Calculate minimal parameter regressor of gravitation load for
% S5RPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% taug_reg [5x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR11_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR11_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t18 = pkin(8) + qJ(3);
t15 = sin(t18);
t16 = cos(t18);
t22 = sin(qJ(5));
t24 = cos(qJ(5));
t36 = -t15 * t24 + t16 * t22;
t23 = sin(qJ(1));
t25 = cos(qJ(1));
t13 = g(1) * t25 + g(2) * t23;
t12 = g(1) * t23 - g(2) * t25;
t31 = t16 * pkin(3) + t15 * qJ(4);
t29 = t15 * t22 + t16 * t24;
t20 = cos(pkin(8));
t28 = t20 * pkin(2) + pkin(1) + t31;
t3 = t36 * t23;
t5 = t36 * t25;
t27 = g(1) * t5 + g(2) * t3 + g(3) * t29;
t4 = t29 * t23;
t6 = t29 * t25;
t26 = g(1) * t6 + g(2) * t4 - g(3) * t36;
t21 = -pkin(6) - qJ(2);
t8 = t12 * t16;
t7 = t12 * t15;
t2 = g(3) * t15 + t13 * t16;
t1 = -g(3) * t16 + t13 * t15;
t9 = [0, t12, t13, t12 * t20, -t12 * sin(pkin(8)), -t13, -g(1) * (-t23 * pkin(1) + t25 * qJ(2)) - g(2) * (t25 * pkin(1) + t23 * qJ(2)), 0, 0, 0, 0, 0, t8, -t7, t8, -t13, t7, (g(1) * t21 - g(2) * t28) * t25 + (g(1) * t28 + g(2) * t21) * t23, 0, 0, 0, 0, 0, g(1) * t4 - g(2) * t6, -g(1) * t3 + g(2) * t5; 0, 0, 0, 0, 0, 0, -t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, 0, -t2, -g(3) * t31 + t13 * (pkin(3) * t15 - qJ(4) * t16), 0, 0, 0, 0, 0, -t27, -t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t26;];
taug_reg = t9;
