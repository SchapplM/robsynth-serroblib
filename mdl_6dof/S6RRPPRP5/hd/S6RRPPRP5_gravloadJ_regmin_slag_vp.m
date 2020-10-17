% Calculate minimal parameter regressor of gravitation load for
% S6RRPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
% 
% Output:
% taug_reg [6x29]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPPRP5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:27:58
% EndTime: 2019-05-06 09:27:59
% DurationCPUTime: 0.34s
% Computational Cost: add. (243->88), mult. (405->108), div. (0->0), fcn. (417->8), ass. (0->55)
t33 = sin(qJ(1));
t35 = cos(qJ(1));
t14 = g(1) * t35 + g(2) * t33;
t32 = sin(qJ(2));
t67 = t14 * t32;
t22 = t32 * qJ(3);
t34 = cos(qJ(2));
t48 = t34 * pkin(2) + t22;
t10 = g(3) * t32 + t14 * t34;
t66 = pkin(2) * t32;
t29 = sin(pkin(9));
t65 = pkin(4) * t29;
t64 = g(1) * t33;
t60 = g(3) * t34;
t28 = pkin(9) + qJ(5);
t20 = sin(t28);
t58 = t33 * t20;
t21 = cos(t28);
t57 = t33 * t21;
t56 = t33 * t29;
t30 = cos(pkin(9));
t55 = t33 * t30;
t54 = t34 * t35;
t53 = t35 * t20;
t52 = t35 * t21;
t51 = t35 * t29;
t50 = t35 * t30;
t47 = qJ(3) * t34;
t46 = t34 * qJ(4);
t45 = t32 * t51;
t44 = g(3) * t48;
t43 = pkin(2) * t54 + t33 * pkin(7) + (pkin(1) + t22) * t35;
t5 = -t32 * t52 + t58;
t7 = t32 * t57 + t53;
t42 = g(1) * t7 + g(2) * t5;
t15 = t33 * t47;
t17 = t35 * t47;
t41 = -g(1) * t17 - g(2) * t15;
t40 = -g(2) * t35 + t64;
t39 = -pkin(1) - t48;
t38 = pkin(5) * t20 - qJ(6) * t21 + t65;
t6 = t32 * t53 + t57;
t8 = -t32 * t58 + t52;
t37 = -g(1) * t6 + g(2) * t8 + t20 * t60;
t1 = g(1) * t5 - g(2) * t7 + t21 * t60;
t31 = -pkin(8) - qJ(4);
t25 = t35 * pkin(7);
t19 = t30 * pkin(4) + pkin(3);
t12 = t40 * t34;
t11 = t40 * t32;
t9 = -t60 + t67;
t4 = t10 * t21;
t3 = t10 * t20;
t2 = -g(1) * t8 - g(2) * t6;
t13 = [0, t40, t14, 0, 0, 0, 0, 0, t12, -t11, -t14, -t12, t11, -g(1) * t25 - g(2) * t43 - t39 * t64, -g(1) * (-t32 * t56 + t50) - g(2) * (t45 + t55) -g(1) * (-t32 * t55 - t51) - g(2) * (t32 * t50 - t56) t12, -g(1) * (t35 * pkin(3) + t25) - g(2) * (t35 * t46 + t43) + (-g(1) * (t39 - t46) - g(2) * pkin(3)) * t33, 0, 0, 0, 0, 0, t2, t42, t2, t12, -t42, -g(1) * (t8 * pkin(5) + t7 * qJ(6) + t35 * t19 + t25) - g(2) * (pkin(4) * t45 + t6 * pkin(5) + t5 * qJ(6) - t31 * t54 + t43) + (-g(1) * (t34 * t31 - t32 * t65 + t39) - g(2) * t19) * t33; 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, -t9, -t10, -g(1) * (-t35 * t66 + t17) - g(2) * (-t33 * t66 + t15) - t44, -t10 * t29, -t10 * t30, t9, -g(3) * (t46 + t48) + t41 + (pkin(2) + qJ(4)) * t67, 0, 0, 0, 0, 0, -t3, -t4, -t3, t9, t4, -t44 + (g(3) * t31 - t14 * t38) * t34 + (-g(3) * t38 + t14 * (pkin(2) - t31)) * t32 + t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t37, t1, 0, t37, -g(1) * (-t5 * pkin(5) + t6 * qJ(6)) - g(2) * (t7 * pkin(5) - t8 * qJ(6)) - (-pkin(5) * t21 - qJ(6) * t20) * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t13;
