% Calculate minimal parameter regressor of gravitation load for
% S5RRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% taug_reg [5x21]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRP6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:58:27
% EndTime: 2019-12-31 19:58:28
% DurationCPUTime: 0.18s
% Computational Cost: add. (115->41), mult. (168->60), div. (0->0), fcn. (163->8), ass. (0->32)
t19 = sin(qJ(1));
t22 = cos(qJ(1));
t6 = g(1) * t22 + g(2) * t19;
t14 = qJ(2) + pkin(8);
t11 = cos(t14);
t20 = cos(qJ(4));
t29 = t22 * t20;
t17 = sin(qJ(4));
t32 = t19 * t17;
t1 = t11 * t32 + t29;
t30 = t22 * t17;
t31 = t19 * t20;
t3 = -t11 * t30 + t31;
t10 = sin(t14);
t34 = g(3) * t10;
t40 = -g(1) * t3 + g(2) * t1 + t17 * t34;
t39 = -g(3) * t11 + t6 * t10;
t16 = -qJ(3) - pkin(6);
t27 = pkin(4) * t17 - t16;
t5 = g(1) * t19 - g(2) * t22;
t15 = -qJ(5) - pkin(7);
t8 = t20 * pkin(4) + pkin(3);
t26 = -t10 * t15 + t11 * t8;
t18 = sin(qJ(2));
t21 = cos(qJ(2));
t23 = -g(3) * t21 + t6 * t18;
t12 = t21 * pkin(2);
t9 = t12 + pkin(1);
t7 = t22 * t9;
t4 = t11 * t29 + t32;
t2 = -t11 * t31 + t30;
t13 = [0, t5, t6, 0, 0, 0, 0, 0, t5 * t21, -t5 * t18, -t6, -g(1) * (-t22 * t16 - t19 * t9) - g(2) * (-t19 * t16 + t7), 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t4, -g(1) * t1 - g(2) * t3, t5 * t10, -g(2) * t7 + (-g(1) * t27 - g(2) * t26) * t22 + (-g(1) * (-t26 - t9) - g(2) * t27) * t19; 0, 0, 0, 0, 0, 0, 0, 0, t23, g(3) * t18 + t6 * t21, 0, t23 * pkin(2), 0, 0, 0, 0, 0, t39 * t20, -t39 * t17, -t6 * t11 - t34, -g(3) * (t12 + t26) + t6 * (pkin(2) * t18 + t10 * t8 + t11 * t15); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, g(1) * t4 - g(2) * t2 + t20 * t34, 0, t40 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39;];
taug_reg = t13;
