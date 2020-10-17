% Calculate minimal parameter regressor of gravitation load for
% S6PRRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% taug_reg [6x29]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRPRR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 05:26:27
% EndTime: 2019-05-05 05:26:28
% DurationCPUTime: 0.38s
% Computational Cost: add. (366->90), mult. (687->159), div. (0->0), fcn. (861->14), ass. (0->49)
t29 = sin(qJ(2));
t31 = cos(qJ(2));
t41 = sin(pkin(11));
t43 = cos(pkin(6));
t36 = t43 * t41;
t42 = cos(pkin(11));
t10 = t42 * t29 + t31 * t36;
t37 = t43 * t42;
t8 = t41 * t29 - t31 * t37;
t56 = -g(1) * t10 - g(2) * t8;
t26 = sin(pkin(6));
t53 = g(3) * t26;
t24 = pkin(12) + qJ(5);
t23 = qJ(6) + t24;
t19 = sin(t23);
t30 = cos(qJ(3));
t52 = t19 * t30;
t20 = cos(t23);
t51 = t20 * t30;
t21 = sin(t24);
t50 = t21 * t30;
t22 = cos(t24);
t49 = t22 * t30;
t25 = sin(pkin(12));
t48 = t25 * t30;
t47 = t26 * t29;
t46 = t26 * t31;
t27 = cos(pkin(12));
t45 = t27 * t30;
t44 = t30 * t31;
t40 = t26 * t42;
t39 = t26 * t41;
t11 = -t29 * t36 + t42 * t31;
t9 = t29 * t37 + t41 * t31;
t38 = -g(1) * t11 - g(2) * t9;
t28 = sin(qJ(3));
t12 = t28 * t47 - t43 * t30;
t4 = t9 * t28 + t30 * t40;
t6 = t11 * t28 - t30 * t39;
t34 = g(1) * t6 + g(2) * t4 + g(3) * t12;
t13 = t43 * t28 + t30 * t47;
t5 = -t28 * t40 + t9 * t30;
t7 = t11 * t30 + t28 * t39;
t33 = g(1) * t7 + g(2) * t5 + g(3) * t13;
t32 = g(3) * t46 + t56;
t3 = t32 * t28;
t2 = -g(1) * (-t10 * t19 - t7 * t20) - g(2) * (-t8 * t19 - t5 * t20) - g(3) * (-t13 * t20 + t19 * t46);
t1 = -g(1) * (t10 * t20 - t7 * t19) - g(2) * (-t5 * t19 + t8 * t20) - g(3) * (-t13 * t19 - t20 * t46);
t14 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t32, g(3) * t47 - t38, 0, 0, 0, 0, 0, -t32 * t30, t3, -g(1) * (-t10 * t45 + t11 * t25) - g(2) * (t9 * t25 - t8 * t45) - (t25 * t29 + t27 * t44) * t53, -g(1) * (t10 * t48 + t11 * t27) - g(2) * (t9 * t27 + t8 * t48) - (-t25 * t44 + t27 * t29) * t53, -t3 (-t29 * t53 + t38) * pkin(8) + (-t31 * t53 - t56) * (pkin(3) * t30 + qJ(4) * t28 + pkin(2)) 0, 0, 0, 0, 0, -g(1) * (-t10 * t49 + t11 * t21) - g(2) * (t9 * t21 - t8 * t49) - (t21 * t29 + t22 * t44) * t53, -g(1) * (t10 * t50 + t11 * t22) - g(2) * (t9 * t22 + t8 * t50) - (-t21 * t44 + t22 * t29) * t53, 0, 0, 0, 0, 0, -g(1) * (-t10 * t51 + t11 * t19) - g(2) * (t9 * t19 - t8 * t51) - (t19 * t29 + t20 * t44) * t53, -g(1) * (t10 * t52 + t11 * t20) - g(2) * (t9 * t20 + t8 * t52) - (-t19 * t44 + t20 * t29) * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t33, t34 * t27, -t34 * t25, -t33, -g(1) * (-t6 * pkin(3) + t7 * qJ(4)) - g(2) * (-t4 * pkin(3) + t5 * qJ(4)) - g(3) * (-t12 * pkin(3) + t13 * qJ(4)) 0, 0, 0, 0, 0, t34 * t22, -t34 * t21, 0, 0, 0, 0, 0, t34 * t20, -t34 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t10 * t22 - t7 * t21) - g(2) * (-t5 * t21 + t8 * t22) - g(3) * (-t13 * t21 - t22 * t46) -g(1) * (-t10 * t21 - t7 * t22) - g(2) * (-t8 * t21 - t5 * t22) - g(3) * (-t13 * t22 + t21 * t46) 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t14;
