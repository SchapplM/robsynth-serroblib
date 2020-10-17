% Calculate minimal parameter regressor of gravitation load for
% S5RRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [5x31]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRR7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR7_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR7_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:22:26
% EndTime: 2019-12-31 22:22:27
% DurationCPUTime: 0.11s
% Computational Cost: add. (174->26), mult. (172->42), div. (0->0), fcn. (180->10), ass. (0->30)
t16 = qJ(2) + qJ(3);
t15 = qJ(4) + t16;
t11 = sin(t15);
t12 = cos(t15);
t19 = sin(qJ(1));
t22 = cos(qJ(1));
t24 = g(1) * t22 + g(2) * t19;
t3 = -g(3) * t12 + t24 * t11;
t30 = g(3) * t11;
t17 = sin(qJ(5));
t28 = t19 * t17;
t20 = cos(qJ(5));
t27 = t19 * t20;
t26 = t22 * t17;
t25 = t22 * t20;
t23 = g(1) * t19 - g(2) * t22;
t21 = cos(qJ(2));
t18 = sin(qJ(2));
t14 = cos(t16);
t13 = sin(t16);
t10 = t12 * t25 + t28;
t9 = -t12 * t26 + t27;
t8 = -t12 * t27 + t26;
t7 = t12 * t28 + t25;
t6 = g(3) * t13 + t24 * t14;
t5 = -g(3) * t14 + t24 * t13;
t4 = t24 * t12 + t30;
t2 = t3 * t20;
t1 = t3 * t17;
t29 = [0, t23, t24, 0, 0, 0, 0, 0, t23 * t21, -t23 * t18, 0, 0, 0, 0, 0, t23 * t14, -t23 * t13, 0, 0, 0, 0, 0, t23 * t12, -t23 * t11, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t10, -g(1) * t7 - g(2) * t9; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t21 + t24 * t18, g(3) * t18 + t24 * t21, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 + g(2) * t7 + t17 * t30, g(1) * t10 - g(2) * t8 + t20 * t30;];
taug_reg = t29;
