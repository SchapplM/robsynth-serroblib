% Calculate minimal parameter regressor of gravitation load for
% S5RRRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [5x23]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRR13_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR13_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 17:32:15
% EndTime: 2024-09-27 17:32:16
% DurationCPUTime: 0.11s
% Computational Cost: add. (388->35), mult. (204->57), div. (0->0), fcn. (216->16), ass. (0->44)
t46 = g(3) * sin(pkin(5));
t34 = cos(pkin(5));
t35 = sin(qJ(4));
t45 = t34 * t35;
t37 = cos(qJ(4));
t44 = t34 * t37;
t32 = qJ(1) + qJ(2);
t31 = qJ(4) + qJ(5);
t25 = pkin(5) + t31;
t43 = sin(t25) / 0.2e1;
t42 = pkin(5) - t31;
t41 = sin(t42);
t15 = t43 - t41 / 0.2e1;
t30 = qJ(3) + t32;
t23 = sin(t30);
t24 = cos(t30);
t28 = cos(t31);
t40 = t24 * t15 + t23 * t28;
t39 = t23 * t15 - t24 * t28;
t38 = cos(qJ(1));
t36 = sin(qJ(1));
t29 = cos(t32);
t27 = sin(t32);
t26 = sin(t31);
t22 = cos(t42);
t21 = cos(t25) / 0.2e1;
t18 = g(1) * t29 + g(2) * t27;
t17 = g(1) * t27 - g(2) * t29;
t16 = t22 / 0.2e1 + t21;
t14 = g(1) * t24 + g(2) * t23;
t13 = g(1) * t23 - g(2) * t24;
t12 = -t23 * t45 + t24 * t37;
t11 = -t23 * t44 - t24 * t35;
t10 = -t23 * t37 - t24 * t45;
t9 = t23 * t35 - t24 * t44;
t8 = -t23 * t16 - t24 * t26;
t7 = -t24 * t16 + t23 * t26;
t6 = -g(1) * t10 - g(2) * t12;
t5 = -g(1) * t9 - g(2) * t11;
t4 = g(1) * t40 + g(2) * t39;
t3 = -g(1) * t7 - g(2) * t8;
t2 = -g(1) * t39 + g(2) * t40 - g(3) * (t21 - t22 / 0.2e1);
t1 = -g(1) * t8 + g(2) * t7 - g(3) * (t43 + t41 / 0.2e1);
t19 = [0, g(1) * t36 - g(2) * t38, g(1) * t38 + g(2) * t36, 0, t17, t18, 0, t13, t14, 0, 0, 0, 0, 0, t6, t5, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, t17, t18, 0, t13, t14, 0, 0, 0, 0, 0, t6, t5, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, t13, t14, 0, 0, 0, 0, 0, t6, t5, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t11 + g(2) * t9 - t37 * t46, g(1) * t12 - g(2) * t10 + t35 * t46, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t19;
