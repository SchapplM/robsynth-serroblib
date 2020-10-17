% Calculate minimal parameter regressor of gravitation load for
% S6PPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
% 
% Output:
% taug_reg [6x20]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PPRPRR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRPRR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_gravloadJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:53:54
% EndTime: 2019-05-04 19:53:55
% DurationCPUTime: 0.30s
% Computational Cost: add. (396->72), mult. (1107->136), div. (0->0), fcn. (1459->16), ass. (0->53)
t33 = sin(pkin(12));
t34 = sin(pkin(11));
t38 = cos(pkin(12));
t39 = cos(pkin(11));
t41 = cos(pkin(6));
t61 = t39 * t41;
t25 = -t34 * t33 + t38 * t61;
t65 = t34 * t41;
t27 = -t39 * t33 - t38 * t65;
t35 = sin(pkin(7));
t40 = cos(pkin(7));
t36 = sin(pkin(6));
t62 = t36 * t40;
t63 = t36 * t39;
t64 = t35 * t36;
t67 = g(1) * (t27 * t40 + t34 * t64) - g(2) * (-t25 * t40 + t35 * t63) + g(3) * (t35 * t41 + t38 * t62);
t66 = t34 * t36;
t42 = sin(qJ(6));
t46 = cos(qJ(5));
t60 = t42 * t46;
t45 = cos(qJ(6));
t59 = t45 * t46;
t32 = sin(pkin(13));
t37 = cos(pkin(13));
t44 = sin(qJ(3));
t47 = cos(qJ(3));
t58 = t47 * t32 + t44 * t37;
t29 = t44 * t32 - t47 * t37;
t16 = -t25 * t35 - t39 * t62;
t17 = -t27 * t35 + t34 * t62;
t24 = -t38 * t64 + t41 * t40;
t43 = sin(qJ(5));
t21 = t58 * t35;
t23 = t58 * t40;
t49 = t41 * t21 + (t23 * t38 - t29 * t33) * t36;
t28 = -t33 * t65 + t39 * t38;
t50 = t21 * t66 + t27 * t23 - t28 * t29;
t26 = t33 * t61 + t34 * t38;
t51 = -t21 * t63 + t25 * t23 - t26 * t29;
t54 = g(1) * (t17 * t46 - t43 * t50) + g(2) * (t16 * t46 - t43 * t51) + g(3) * (t24 * t46 - t43 * t49);
t20 = t29 * t35;
t22 = t29 * t40;
t14 = -t41 * t20 + (-t22 * t38 - t33 * t58) * t36;
t6 = t20 * t63 - t25 * t22 - t26 * t58;
t9 = -t20 * t66 - t27 * t22 - t28 * t58;
t53 = g(1) * t9 + g(2) * t6 + g(3) * t14;
t52 = g(3) * t33 * t36 + g(1) * t28 + g(2) * t26;
t48 = t52 * t44 - t67 * t47;
t19 = -g(3) * t41 + (-g(1) * t34 + g(2) * t39) * t36;
t12 = t24 * t43 + t46 * t49;
t4 = t17 * t43 + t46 * t50;
t2 = t16 * t43 + t46 * t51;
t1 = [-g(3), -g(3), 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, t19, 0, 0, 0, t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, t48, t67 * t44 + t52 * t47, t48 * pkin(3), 0, 0, 0, 0, 0, -t53 * t46, t53 * t43, 0, 0, 0, 0, 0, -g(1) * (t42 * t50 + t9 * t59) - g(2) * (t42 * t51 + t6 * t59) - g(3) * (t14 * t59 + t42 * t49) -g(1) * (t45 * t50 - t9 * t60) - g(2) * (t45 * t51 - t6 * t60) - g(3) * (-t14 * t60 + t45 * t49); 0, 0, 0, 0, 0, -g(1) * t17 - g(2) * t16 - g(3) * t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, g(1) * t4 + g(2) * t2 + g(3) * t12, 0, 0, 0, 0, 0, -t54 * t45, t54 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t4 * t42 - t9 * t45) - g(2) * (-t2 * t42 - t6 * t45) - g(3) * (-t12 * t42 - t14 * t45) -g(1) * (-t4 * t45 + t9 * t42) - g(2) * (-t2 * t45 + t6 * t42) - g(3) * (-t12 * t45 + t14 * t42);];
taug_reg  = t1;
