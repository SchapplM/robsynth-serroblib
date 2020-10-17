% Calculate minimal parameter regressor of gravitation load for
% S6RRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% 
% Output:
% taug_reg [6x28]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP1_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 07:27:09
% EndTime: 2019-05-07 07:27:10
% DurationCPUTime: 0.25s
% Computational Cost: add. (260->52), mult. (256->72), div. (0->0), fcn. (245->10), ass. (0->45)
t26 = qJ(2) + qJ(3);
t20 = pkin(10) + t26;
t16 = sin(t20);
t17 = cos(t20);
t31 = cos(qJ(5));
t19 = t31 * pkin(5) + pkin(4);
t27 = -qJ(6) - pkin(9);
t38 = -t16 * t27 + t17 * t19;
t30 = sin(qJ(1));
t33 = cos(qJ(1));
t15 = g(1) * t33 + g(2) * t30;
t28 = sin(qJ(5));
t49 = g(3) * t16;
t43 = t33 * t31;
t46 = t30 * t28;
t6 = t17 * t46 + t43;
t44 = t33 * t28;
t45 = t30 * t31;
t8 = -t17 * t44 + t45;
t56 = -g(1) * t8 + g(2) * t6 + t28 * t49;
t55 = -g(3) * t17 + t15 * t16;
t21 = sin(t26);
t52 = pkin(3) * t21;
t22 = cos(t26);
t18 = pkin(3) * t22;
t32 = cos(qJ(2));
t23 = t32 * pkin(2);
t42 = t18 + t23;
t25 = -qJ(4) - pkin(8) - pkin(7);
t40 = pkin(5) * t28 - t25;
t39 = t18 + t38;
t14 = g(1) * t30 - g(2) * t33;
t37 = -t16 * t19 - t17 * t27;
t4 = -g(3) * t22 + t15 * t21;
t29 = sin(qJ(2));
t13 = -t29 * pkin(2) - t52;
t12 = pkin(1) + t42;
t10 = t33 * t12;
t9 = t17 * t43 + t46;
t7 = -t17 * t45 + t44;
t5 = g(3) * t21 + t15 * t22;
t3 = -t15 * t17 - t49;
t2 = t55 * t31;
t1 = t55 * t28;
t11 = [0, t14, t15, 0, 0, 0, 0, 0, t14 * t32, -t14 * t29, 0, 0, 0, 0, 0, t14 * t22, -t14 * t21, -t15, -g(1) * (-t30 * t12 - t33 * t25) - g(2) * (-t30 * t25 + t10) 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t9, -g(1) * t6 - g(2) * t8, t14 * t16, -g(2) * t10 + (-g(1) * t40 - g(2) * t38) * t33 + (-g(1) * (-t12 - t38) - g(2) * t40) * t30; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t32 + t15 * t29, g(3) * t29 + t15 * t32, 0, 0, 0, 0, 0, t4, t5, 0, -g(3) * t42 - t15 * t13, 0, 0, 0, 0, 0, t2, -t1, t3, -g(3) * (t23 + t39) + t15 * (-t13 - t37); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t5, 0, t4 * pkin(3), 0, 0, 0, 0, 0, t2, -t1, t3, -g(3) * t39 + t15 * (-t37 + t52); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, g(1) * t9 - g(2) * t7 + t31 * t49, 0, t56 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55;];
taug_reg  = t11;
