% Calculate minimal parameter regressor of gravitation load for
% S6PRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
% 
% Output:
% taug_reg [6x24]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPPRR3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:01:35
% EndTime: 2019-05-04 22:01:36
% DurationCPUTime: 0.25s
% Computational Cost: add. (218->70), mult. (567->118), div. (0->0), fcn. (725->12), ass. (0->41)
t36 = sin(pkin(6));
t40 = sin(qJ(5));
t61 = t36 * t40;
t41 = sin(qJ(2));
t60 = t36 * t41;
t43 = cos(qJ(5));
t59 = t36 * t43;
t44 = cos(qJ(2));
t58 = t36 * t44;
t39 = sin(qJ(6));
t57 = t39 * t43;
t42 = cos(qJ(6));
t56 = t42 * t43;
t55 = pkin(2) * t58 + qJ(3) * t60;
t54 = cos(pkin(6));
t35 = sin(pkin(10));
t38 = cos(pkin(10));
t50 = t44 * t54;
t25 = t35 * t41 - t38 * t50;
t51 = t41 * t54;
t26 = t35 * t44 + t38 * t51;
t34 = sin(pkin(11));
t37 = cos(pkin(11));
t53 = -t25 * t37 + t26 * t34;
t10 = t25 * t34 + t26 * t37;
t27 = t35 * t50 + t38 * t41;
t28 = -t35 * t51 + t38 * t44;
t52 = -t27 * t37 + t28 * t34;
t14 = t27 * t34 + t28 * t37;
t49 = -t25 * pkin(2) + t26 * qJ(3);
t48 = -t27 * pkin(2) + t28 * qJ(3);
t22 = -t34 * t58 + t37 * t60;
t47 = g(1) * (-t14 * t40 - t35 * t59) + g(2) * (-t10 * t40 + t38 * t59) + g(3) * (-t22 * t40 - t54 * t43);
t21 = (t34 * t41 + t37 * t44) * t36;
t46 = g(1) * t52 + g(2) * t53 + g(3) * t21;
t5 = -g(1) * t27 - g(2) * t25 + g(3) * t58;
t45 = g(1) * t28 + g(2) * t26 + g(3) * t60;
t16 = t22 * t43 - t54 * t40;
t4 = t14 * t43 - t35 * t61;
t2 = t10 * t43 + t38 * t61;
t1 = [-g(3), 0, 0, 0, 0, 0, -g(3), 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t5, t45, -t5, -t45, -g(1) * t48 - g(2) * t49 - g(3) * t55, -t46, -g(1) * t14 - g(2) * t10 - g(3) * t22, -g(1) * (-t27 * pkin(3) + t48) - g(2) * (-t25 * pkin(3) + t49) - g(3) * (pkin(3) * t58 + t55) 0, 0, 0, 0, 0, -t46 * t43, t46 * t40, 0, 0, 0, 0, 0, -g(1) * (-t14 * t39 + t52 * t56) - g(2) * (-t10 * t39 + t53 * t56) - g(3) * (t21 * t56 - t22 * t39) -g(1) * (-t14 * t42 - t52 * t57) - g(2) * (-t10 * t42 - t53 * t57) - g(3) * (-t21 * t57 - t22 * t42); 0, 0, 0, 0, 0, 0, t5, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t54 + (g(1) * t35 - g(2) * t38) * t36, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, g(1) * t4 + g(2) * t2 + g(3) * t16, 0, 0, 0, 0, 0, -t47 * t42, t47 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t4 * t39 + t42 * t52) - g(2) * (-t2 * t39 + t42 * t53) - g(3) * (-t16 * t39 + t21 * t42) -g(1) * (-t39 * t52 - t4 * t42) - g(2) * (-t2 * t42 - t39 * t53) - g(3) * (-t16 * t42 - t21 * t39);];
taug_reg  = t1;
