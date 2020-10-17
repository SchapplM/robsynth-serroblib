% Calculate minimal parameter regressor of gravitation load for
% S6RRPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
% 
% Output:
% taug_reg [6x29]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPPPR3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 08:35:27
% EndTime: 2019-05-06 08:35:28
% DurationCPUTime: 0.29s
% Computational Cost: add. (147->71), mult. (300->95), div. (0->0), fcn. (297->8), ass. (0->53)
t28 = sin(qJ(2));
t18 = t28 * qJ(3);
t62 = pkin(1) + t18;
t31 = cos(qJ(1));
t29 = sin(qJ(1));
t57 = g(2) * t29;
t10 = g(1) * t31 + t57;
t61 = t10 * t28;
t30 = cos(qJ(2));
t6 = g(3) * t28 + t10 * t30;
t60 = pkin(2) + pkin(3);
t59 = g(1) * t29;
t55 = g(3) * t30;
t21 = t30 * pkin(2);
t20 = t30 * pkin(3);
t54 = g(2) * qJ(4);
t53 = t28 * t31;
t25 = pkin(9) + qJ(6);
t16 = sin(t25);
t52 = t29 * t16;
t17 = cos(t25);
t51 = t29 * t17;
t26 = sin(pkin(9));
t50 = t29 * t26;
t27 = cos(pkin(9));
t49 = t29 * t27;
t48 = t30 * t31;
t47 = t31 * t16;
t46 = t31 * t17;
t45 = t31 * t26;
t44 = t31 * t27;
t43 = t21 + t18;
t42 = qJ(3) * t30;
t41 = t30 * qJ(5);
t40 = qJ(5) + t60;
t39 = t20 + t43;
t37 = pkin(2) * t48 + t29 * pkin(7) + t31 * t62;
t36 = g(1) * t40;
t35 = pkin(3) * t48 + t37;
t9 = -g(2) * t31 + t59;
t22 = t31 * pkin(7);
t34 = g(1) * (-t31 * qJ(4) + t22);
t33 = -t62 - t21;
t13 = t31 * t42;
t11 = t29 * t42;
t8 = t9 * t30;
t7 = t9 * t28;
t5 = -t55 + t61;
t4 = t28 * t46 - t52;
t3 = -t28 * t47 - t51;
t2 = -t28 * t51 - t47;
t1 = t28 * t52 - t46;
t12 = [0, t9, t10, 0, 0, 0, 0, 0, t8, -t7, t8, -t10, t7, -g(1) * t22 - g(2) * t37 - t33 * t59, t7, -t8, t10, -t34 - g(2) * t35 + (-g(1) * (t33 - t20) + t54) * t29, -g(1) * (-t28 * t49 - t45) - g(2) * (t28 * t44 - t50) -g(1) * (t28 * t50 - t44) - g(2) * (-t28 * t45 - t49) t8, -t34 - g(2) * (pkin(4) * t53 + t31 * t41 + t35) + (-g(1) * (-t28 * pkin(4) - t62) + t54 + t30 * t36) * t29, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t4, -g(1) * t1 - g(2) * t3; 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, t5, 0, -t6, -g(1) * (-pkin(2) * t53 + t13) - g(2) * (-pkin(2) * t28 * t29 + t11) - g(3) * t43, -t6, -t5, 0, -g(1) * t13 - g(2) * t11 - g(3) * t39 + t60 * t61, -t6 * t27, t6 * t26, t5, -g(1) * (pkin(4) * t48 + t13) - g(2) * (t29 * t30 * pkin(4) + t11) - g(3) * (t39 + t41) + (-g(3) * pkin(4) + t31 * t36 + t40 * t57) * t28, 0, 0, 0, 0, 0, -t6 * t17, t6 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, -t5, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 + g(2) * t1 - t16 * t55, g(1) * t4 - g(2) * t2 - t17 * t55;];
taug_reg  = t12;
