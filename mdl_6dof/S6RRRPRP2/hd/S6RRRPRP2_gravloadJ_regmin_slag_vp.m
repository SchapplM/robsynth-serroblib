% Calculate minimal parameter regressor of gravitation load for
% S6RRRPRP2
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
% taug_reg [6x30]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRP2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP2_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP2_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 07:32:55
% EndTime: 2019-05-07 07:32:56
% DurationCPUTime: 0.28s
% Computational Cost: add. (373->70), mult. (373->95), div. (0->0), fcn. (380->10), ass. (0->52)
t36 = sin(qJ(5));
t39 = cos(qJ(5));
t63 = pkin(5) * t39 + qJ(6) * t36;
t38 = sin(qJ(1));
t41 = cos(qJ(1));
t21 = g(1) * t41 + g(2) * t38;
t35 = qJ(2) + qJ(3);
t29 = pkin(10) + t35;
t24 = sin(t29);
t44 = t21 * t24;
t22 = t24 * pkin(9);
t25 = cos(t29);
t62 = -pkin(4) * t25 - t22;
t30 = sin(t35);
t61 = pkin(3) * t30;
t57 = g(3) * t24;
t56 = g(3) * t36;
t55 = t25 * t41;
t54 = t38 * t36;
t53 = t38 * t39;
t34 = -qJ(4) - pkin(8) - pkin(7);
t52 = t41 * t34;
t51 = t41 * t36;
t50 = t41 * t39;
t31 = cos(t35);
t26 = pkin(3) * t31;
t40 = cos(qJ(2));
t32 = t40 * pkin(2);
t49 = t26 + t32;
t10 = t25 * t51 - t53;
t8 = t25 * t54 + t50;
t47 = g(1) * t8 - g(2) * t10;
t46 = t25 * t63 + t26 - t62;
t20 = g(1) * t38 - g(2) * t41;
t1 = g(1) * t10 + g(2) * t8 + t24 * t56;
t11 = t25 * t50 + t54;
t9 = t25 * t53 - t51;
t43 = g(1) * t11 + g(2) * t9 + t39 * t57;
t6 = -g(3) * t31 + t21 * t30;
t42 = (pkin(4) + t63) * t44;
t37 = sin(qJ(2));
t18 = pkin(9) * t55;
t16 = t38 * t25 * pkin(9);
t14 = -pkin(2) * t37 - t61;
t13 = pkin(1) + t49;
t12 = t41 * t13;
t7 = g(3) * t30 + t21 * t31;
t5 = -t21 * t25 - t57;
t4 = (-g(3) * t25 + t44) * t39;
t3 = -t25 * t56 + t36 * t44;
t2 = g(1) * t9 - g(2) * t11;
t15 = [0, t20, t21, 0, 0, 0, 0, 0, t20 * t40, -t20 * t37, 0, 0, 0, 0, 0, t20 * t31, -t20 * t30, -t21, -g(1) * (-t13 * t38 - t52) - g(2) * (-t38 * t34 + t12) 0, 0, 0, 0, 0, t2, -t47, t2, t20 * t24, t47, -g(1) * (-t9 * pkin(5) - t8 * qJ(6) - t52) - g(2) * (pkin(4) * t55 + t11 * pkin(5) + t10 * qJ(6) + t22 * t41 + t12) + (-g(1) * (-t13 + t62) + g(2) * t34) * t38; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t40 + t21 * t37, g(3) * t37 + t21 * t40, 0, 0, 0, 0, 0, t6, t7, 0, -g(3) * t49 - t14 * t21, 0, 0, 0, 0, 0, t4, -t3, t4, t5, t3, -g(1) * (t41 * t14 + t18) - g(2) * (t38 * t14 + t16) - g(3) * (t32 + t46) + t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t7, 0, t6 * pkin(3), 0, 0, 0, 0, 0, t4, -t3, t4, t5, t3, -g(1) * (-t41 * t61 + t18) - g(2) * (-t38 * t61 + t16) - g(3) * t46 + t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t43, t1, 0, -t43, -g(1) * (-pkin(5) * t10 + qJ(6) * t11) - g(2) * (-pkin(5) * t8 + qJ(6) * t9) - (-pkin(5) * t36 + qJ(6) * t39) * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t15;
