% Calculate minimal parameter regressor of gravitation load for
% S5PRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% taug_reg [5x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRRR11_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR11_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 21:45:57
% EndTime: 2024-09-27 21:45:57
% DurationCPUTime: 0.19s
% Computational Cost: add. (423->48), mult. (210->73), div. (0->0), fcn. (208->18), ass. (0->48)
t52 = g(3) * sin(pkin(5));
t37 = cos(pkin(5));
t38 = sin(qJ(3));
t51 = t37 * t38;
t39 = cos(qJ(3));
t50 = t37 * t39;
t35 = qJ(3) + qJ(4);
t30 = pkin(5) + t35;
t25 = qJ(5) + t30;
t49 = sin(t25) / 0.2e1;
t48 = sin(t30) / 0.2e1;
t47 = pkin(5) - t35;
t46 = -qJ(5) + t47;
t45 = sin(t47);
t40 = sin(t46);
t13 = t49 - t40 / 0.2e1;
t33 = qJ(5) + t35;
t27 = cos(t33);
t34 = pkin(10) + qJ(2);
t28 = sin(t34);
t29 = cos(t34);
t44 = t29 * t13 + t28 * t27;
t43 = t28 * t13 - t29 * t27;
t15 = t48 - t45 / 0.2e1;
t32 = cos(t35);
t42 = t29 * t15 + t28 * t32;
t41 = t28 * t15 - t29 * t32;
t31 = sin(t35);
t26 = sin(t33);
t24 = cos(t47);
t23 = cos(t46);
t22 = cos(t30) / 0.2e1;
t19 = cos(t25) / 0.2e1;
t16 = t24 / 0.2e1 + t22;
t14 = t23 / 0.2e1 + t19;
t12 = -t28 * t51 + t29 * t39;
t11 = -t28 * t50 - t29 * t38;
t10 = -t28 * t39 - t29 * t51;
t9 = t28 * t38 - t29 * t50;
t8 = -t28 * t16 - t29 * t31;
t7 = -t29 * t16 + t28 * t31;
t6 = -t28 * t14 - t29 * t26;
t5 = -t29 * t14 + t28 * t26;
t4 = -g(1) * t41 + g(2) * t42 - g(3) * (t22 - t24 / 0.2e1);
t3 = -g(1) * t8 + g(2) * t7 - g(3) * (t48 + t45 / 0.2e1);
t2 = -g(1) * t43 + g(2) * t44 - g(3) * (t19 - t23 / 0.2e1);
t1 = -g(1) * t6 + g(2) * t5 - g(3) * (t49 + t40 / 0.2e1);
t17 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, g(1) * t28 - g(2) * t29, g(1) * t29 + g(2) * t28, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, -g(1) * t9 - g(2) * t11, 0, 0, 0, 0, 0, g(1) * t42 + g(2) * t41, -g(1) * t7 - g(2) * t8, 0, 0, 0, 0, 0, g(1) * t44 + g(2) * t43, -g(1) * t5 - g(2) * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t11 + g(2) * t9 - t39 * t52, g(1) * t12 - g(2) * t10 + t38 * t52, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t17;
