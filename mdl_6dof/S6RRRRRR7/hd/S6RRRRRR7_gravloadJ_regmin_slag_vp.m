% Calculate minimal parameter regressor of gravitation load for
% S6RRRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% taug_reg [6x38]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 04:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRR7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR7_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR7_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 12:33:23
% EndTime: 2019-05-08 12:33:26
% DurationCPUTime: 0.59s
% Computational Cost: add. (516->98), mult. (928->187), div. (0->0), fcn. (1198->14), ass. (0->58)
t33 = sin(qJ(2));
t34 = sin(qJ(1));
t37 = cos(qJ(2));
t44 = cos(pkin(6));
t61 = cos(qJ(1));
t41 = t44 * t61;
t19 = t33 * t41 + t34 * t37;
t32 = sin(qJ(3));
t36 = cos(qJ(3));
t30 = sin(pkin(6));
t43 = t30 * t61;
t12 = t19 * t36 - t32 * t43;
t18 = t34 * t33 - t37 * t41;
t29 = qJ(4) + qJ(5);
t28 = qJ(6) + t29;
t24 = sin(t28);
t25 = cos(t28);
t68 = t12 * t24 - t18 * t25;
t67 = t12 * t25 + t18 * t24;
t26 = sin(t29);
t27 = cos(t29);
t66 = t12 * t26 - t18 * t27;
t65 = t12 * t27 + t18 * t26;
t31 = sin(qJ(4));
t35 = cos(qJ(4));
t64 = t12 * t31 - t18 * t35;
t63 = t12 * t35 + t18 * t31;
t62 = g(3) * t30;
t54 = t24 * t36;
t53 = t25 * t36;
t52 = t26 * t36;
t51 = t27 * t36;
t50 = t30 * t33;
t49 = t30 * t36;
t48 = t30 * t37;
t47 = t31 * t36;
t46 = t35 * t36;
t45 = t36 * t37;
t42 = t34 * t44;
t21 = -t33 * t42 + t61 * t37;
t14 = -t21 * t32 + t34 * t49;
t39 = t19 * t32 + t36 * t43;
t40 = g(1) * t14 - g(2) * t39 + g(3) * (-t32 * t50 + t44 * t36);
t20 = t61 * t33 + t37 * t42;
t38 = -g(1) * t20 - g(2) * t18 + g(3) * t48;
t17 = t44 * t32 + t33 * t49;
t15 = t34 * t30 * t32 + t21 * t36;
t10 = t15 * t35 + t20 * t31;
t9 = -t15 * t31 + t20 * t35;
t8 = t15 * t27 + t20 * t26;
t7 = -t15 * t26 + t20 * t27;
t6 = t15 * t25 + t20 * t24;
t5 = -t15 * t24 + t20 * t25;
t4 = g(1) * t8 + g(2) * t65 - g(3) * (-t17 * t27 + t26 * t48);
t3 = -g(1) * t7 + g(2) * t66 - g(3) * (-t17 * t26 - t27 * t48);
t2 = g(1) * t6 + g(2) * t67 - g(3) * (-t17 * t25 + t24 * t48);
t1 = -g(1) * t5 + g(2) * t68 - g(3) * (-t17 * t24 - t25 * t48);
t11 = [0, g(1) * t34 - g(2) * t61, g(1) * t61 + g(2) * t34, 0, 0, 0, 0, 0, g(1) * t19 - g(2) * t21, -g(1) * t18 + g(2) * t20, 0, 0, 0, 0, 0, g(1) * t12 - g(2) * t15, -g(1) * t39 - g(2) * t14, 0, 0, 0, 0, 0, g(1) * t63 - g(2) * t10, -g(1) * t64 - g(2) * t9, 0, 0, 0, 0, 0, g(1) * t65 - g(2) * t8, -g(1) * t66 - g(2) * t7, 0, 0, 0, 0, 0, g(1) * t67 - g(2) * t6, -g(1) * t68 - g(2) * t5; 0, 0, 0, 0, 0, 0, 0, 0, -t38, g(1) * t21 + g(2) * t19 + g(3) * t50, 0, 0, 0, 0, 0, -t38 * t36, t38 * t32, 0, 0, 0, 0, 0, -g(1) * (-t20 * t46 + t21 * t31) - g(2) * (-t18 * t46 + t19 * t31) - (t31 * t33 + t35 * t45) * t62, -g(1) * (t20 * t47 + t21 * t35) - g(2) * (t18 * t47 + t19 * t35) - (-t31 * t45 + t33 * t35) * t62, 0, 0, 0, 0, 0, -g(1) * (-t20 * t51 + t21 * t26) - g(2) * (-t18 * t51 + t19 * t26) - (t26 * t33 + t27 * t45) * t62, -g(1) * (t20 * t52 + t21 * t27) - g(2) * (t18 * t52 + t19 * t27) - (-t26 * t45 + t27 * t33) * t62, 0, 0, 0, 0, 0, -g(1) * (-t20 * t53 + t21 * t24) - g(2) * (-t18 * t53 + t19 * t24) - (t24 * t33 + t25 * t45) * t62, -g(1) * (t20 * t54 + t21 * t25) - g(2) * (t18 * t54 + t19 * t25) - (-t24 * t45 + t25 * t33) * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, g(1) * t15 + g(2) * t12 + g(3) * t17, 0, 0, 0, 0, 0, -t40 * t35, t40 * t31, 0, 0, 0, 0, 0, -t40 * t27, t40 * t26, 0, 0, 0, 0, 0, -t40 * t25, t40 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 + g(2) * t64 - g(3) * (-t17 * t31 - t35 * t48) g(1) * t10 + g(2) * t63 - g(3) * (-t17 * t35 + t31 * t48) 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t11;
