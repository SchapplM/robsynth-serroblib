% Calculate minimal parameter regressor of gravitation load for
% S5RRPPR5
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
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPPR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t17 = qJ(2) + pkin(8);
t13 = sin(t17);
t14 = cos(t17);
t19 = sin(qJ(5));
t22 = cos(qJ(5));
t35 = -t13 * t22 + t14 * t19;
t21 = sin(qJ(1));
t24 = cos(qJ(1));
t10 = g(1) * t24 + g(2) * t21;
t9 = g(1) * t21 - g(2) * t24;
t30 = t14 * pkin(3) + t13 * qJ(4);
t29 = t13 * t19 + t14 * t22;
t2 = t35 * t21;
t4 = t35 * t24;
t28 = g(1) * t4 + g(2) * t2 + g(3) * t29;
t3 = t29 * t21;
t5 = t29 * t24;
t27 = g(1) * t5 + g(2) * t3 - g(3) * t35;
t20 = sin(qJ(2));
t23 = cos(qJ(2));
t25 = -g(3) * t23 + t10 * t20;
t18 = -qJ(3) - pkin(6);
t15 = t23 * pkin(2);
t12 = t15 + pkin(1);
t11 = t24 * t12;
t1 = -g(3) * t14 + t10 * t13;
t6 = [0, t9, t10, 0, 0, 0, 0, 0, t9 * t23, -t9 * t20, -t10, -g(1) * (-t21 * t12 - t24 * t18) - g(2) * (-t21 * t18 + t11), t9 * t14, -t10, t9 * t13, -g(2) * t11 + (g(1) * t18 - g(2) * t30) * t24 + (-g(1) * (-t12 - t30) + g(2) * t18) * t21, 0, 0, 0, 0, 0, g(1) * t3 - g(2) * t5, -g(1) * t2 + g(2) * t4; 0, 0, 0, 0, 0, 0, 0, 0, t25, g(3) * t20 + t10 * t23, 0, t25 * pkin(2), t1, 0, -g(3) * t13 - t10 * t14, -g(3) * (t15 + t30) + t10 * (pkin(2) * t20 + pkin(3) * t13 - qJ(4) * t14), 0, 0, 0, 0, 0, -t28, -t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t27;];
taug_reg = t6;
