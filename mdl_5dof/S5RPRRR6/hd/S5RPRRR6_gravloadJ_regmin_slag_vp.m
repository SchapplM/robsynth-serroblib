% Calculate minimal parameter regressor of gravitation load for
% S5RPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% taug_reg [5x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRR6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:01:40
% EndTime: 2019-12-31 19:01:40
% DurationCPUTime: 0.11s
% Computational Cost: add. (124->25), mult. (122->43), div. (0->0), fcn. (128->10), ass. (0->27)
t14 = qJ(3) + qJ(4);
t11 = sin(t14);
t12 = cos(t14);
t13 = qJ(1) + pkin(9);
t10 = cos(t13);
t9 = sin(t13);
t22 = g(1) * t10 + g(2) * t9;
t3 = -g(3) * t12 + t22 * t11;
t27 = g(3) * t11;
t15 = sin(qJ(5));
t25 = t12 * t15;
t18 = cos(qJ(5));
t24 = t12 * t18;
t23 = g(1) * t9 - g(2) * t10;
t17 = sin(qJ(1));
t20 = cos(qJ(1));
t21 = g(1) * t17 - g(2) * t20;
t19 = cos(qJ(3));
t16 = sin(qJ(3));
t8 = t10 * t24 + t9 * t15;
t7 = -t10 * t25 + t9 * t18;
t6 = t10 * t15 - t9 * t24;
t5 = t10 * t18 + t9 * t25;
t4 = t22 * t12 + t27;
t2 = t3 * t18;
t1 = t3 * t15;
t26 = [0, t21, g(1) * t20 + g(2) * t17, t21 * pkin(1), 0, 0, 0, 0, 0, t23 * t19, -t23 * t16, 0, 0, 0, 0, 0, t23 * t12, -t23 * t11, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t8, -g(1) * t5 - g(2) * t7; 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t19 + t22 * t16, g(3) * t16 + t22 * t19, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 + g(2) * t5 + t15 * t27, g(1) * t8 - g(2) * t6 + t18 * t27;];
taug_reg = t26;
