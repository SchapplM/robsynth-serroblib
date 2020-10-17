% Calculate minimal parameter regressor of gravitation load for
% S6PRRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
% 
% Output:
% taug_reg [6x29]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRPR3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 07:30:43
% EndTime: 2019-05-05 07:30:44
% DurationCPUTime: 0.29s
% Computational Cost: add. (386->77), mult. (655->119), div. (0->0), fcn. (803->12), ass. (0->46)
t31 = sin(pkin(11));
t35 = sin(qJ(2));
t38 = cos(qJ(2));
t49 = cos(pkin(11));
t50 = cos(pkin(6));
t45 = t50 * t49;
t18 = t31 * t35 - t38 * t45;
t48 = t31 * t50;
t20 = t49 * t35 + t38 * t48;
t32 = sin(pkin(6));
t58 = g(3) * t32;
t43 = -g(1) * t20 - g(2) * t18 + t38 * t58;
t30 = qJ(3) + qJ(4);
t28 = sin(t30);
t33 = sin(qJ(6));
t57 = t28 * t33;
t36 = cos(qJ(6));
t56 = t28 * t36;
t55 = t31 * t32;
t54 = t32 * t35;
t37 = cos(qJ(3));
t53 = t32 * t37;
t52 = t33 * t38;
t51 = t36 * t38;
t47 = t32 * t49;
t19 = t31 * t38 + t35 * t45;
t21 = -t35 * t48 + t49 * t38;
t46 = g(1) * t21 + g(2) * t19;
t29 = cos(t30);
t11 = t19 * t28 + t29 * t47;
t13 = t21 * t28 - t29 * t55;
t16 = t28 * t54 - t50 * t29;
t4 = g(1) * t13 + g(2) * t11 + g(3) * t16;
t12 = t19 * t29 - t28 * t47;
t14 = t21 * t29 + t28 * t55;
t17 = t50 * t28 + t29 * t54;
t6 = g(1) * t14 + g(2) * t12 + g(3) * t17;
t42 = g(3) * t54 + t46;
t41 = -g(1) * (-t13 * pkin(4) + t14 * qJ(5)) - g(2) * (-t11 * pkin(4) + t12 * qJ(5)) - g(3) * (-t16 * pkin(4) + t17 * qJ(5));
t34 = sin(qJ(3));
t40 = -g(1) * (-t21 * t34 + t31 * t53) - g(2) * (-t19 * t34 - t37 * t47) - g(3) * (-t34 * t54 + t50 * t37);
t8 = t43 * t29;
t7 = t43 * t28;
t2 = t6 * t36;
t1 = t6 * t33;
t3 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, -t43, t42, 0, 0, 0, 0, 0, -t43 * t37, t43 * t34, 0, 0, 0, 0, 0, -t8, t7, -t42, t8, -t7 (t35 * t58 + t46) * (-pkin(9) - pkin(8)) - t43 * (t37 * pkin(3) + pkin(4) * t29 + qJ(5) * t28 + pkin(2)) 0, 0, 0, 0, 0, -g(1) * (-t20 * t57 + t21 * t36) - g(2) * (-t18 * t57 + t19 * t36) - (t28 * t52 + t35 * t36) * t58, -g(1) * (-t20 * t56 - t21 * t33) - g(2) * (-t18 * t56 - t19 * t33) - (t28 * t51 - t33 * t35) * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -g(1) * (-t21 * t37 - t34 * t55) - g(2) * (-t19 * t37 + t34 * t47) - g(3) * (-t50 * t34 - t35 * t53) 0, 0, 0, 0, 0, t4, t6, 0, -t4, -t6, t40 * pkin(3) + t41, 0, 0, 0, 0, 0, -t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t6, 0, -t4, -t6, t41, 0, 0, 0, 0, 0, -t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t13 * t36 - t20 * t33) - g(2) * (t11 * t36 - t18 * t33) - g(3) * (t16 * t36 + t32 * t52) -g(1) * (-t13 * t33 - t20 * t36) - g(2) * (-t11 * t33 - t18 * t36) - g(3) * (-t16 * t33 + t32 * t51);];
taug_reg  = t3;
