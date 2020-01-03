% Calculate minimal parameter regressor of gravitation load for
% S5RPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% taug_reg [5x18]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t19 = sin(pkin(9));
t30 = g(1) * t19;
t20 = cos(pkin(9));
t21 = sin(qJ(5));
t29 = t20 * t21;
t23 = cos(qJ(5));
t28 = t20 * t23;
t18 = qJ(1) + pkin(8);
t17 = qJ(3) + t18;
t15 = sin(t17);
t16 = cos(t17);
t27 = t16 * pkin(3) + t15 * qJ(4);
t26 = t15 * pkin(3) - t16 * qJ(4);
t10 = g(2) * t16 + g(3) * t15;
t22 = sin(qJ(1));
t24 = cos(qJ(1));
t25 = -g(2) * t24 - g(3) * t22;
t9 = g(2) * t15 - g(3) * t16;
t8 = t10 * t20;
t7 = t10 * t19;
t6 = t15 * t21 + t16 * t28;
t5 = -t15 * t23 + t16 * t29;
t4 = t15 * t28 - t16 * t21;
t3 = -t15 * t29 - t16 * t23;
t2 = -g(2) * t6 - g(3) * t4;
t1 = g(2) * t5 - g(3) * t3;
t11 = [0, t25, g(2) * t22 - g(3) * t24, t25 * pkin(1), 0, -t10, t9, -t8, t7, -t9, -g(2) * (pkin(2) * cos(t18) + t24 * pkin(1) + t27) - g(3) * (pkin(2) * sin(t18) + t22 * pkin(1) + t26), 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t10, t9, -t8, t7, -t9, -g(2) * t27 - g(3) * t26, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2) * t3 - g(3) * t5 + t21 * t30, g(2) * t4 - g(3) * t6 + t23 * t30;];
taug_reg = t11;
