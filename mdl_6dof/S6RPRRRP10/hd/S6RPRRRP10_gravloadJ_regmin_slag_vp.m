% Calculate minimal parameter regressor of gravitation load for
% S6RPRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% 
% Output:
% taug_reg [6x31]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRP10_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP10_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:58:03
% EndTime: 2019-05-06 01:58:03
% DurationCPUTime: 0.29s
% Computational Cost: add. (262->67), mult. (379->91), div. (0->0), fcn. (412->8), ass. (0->50)
t33 = sin(qJ(1));
t36 = cos(qJ(1));
t67 = -g(1) * t33 + g(2) * t36;
t32 = sin(qJ(3));
t34 = cos(qJ(4));
t47 = t36 * t34;
t31 = sin(qJ(4));
t53 = t33 * t31;
t14 = -t32 * t53 + t47;
t48 = t36 * t31;
t52 = t33 * t34;
t16 = t32 * t48 + t52;
t35 = cos(qJ(3));
t59 = g(3) * t35;
t66 = -g(1) * t14 - g(2) * t16 + t31 * t59;
t40 = -g(3) * t32 - t67 * t35;
t30 = qJ(4) + qJ(5);
t24 = sin(t30);
t58 = t24 * t35;
t25 = cos(t30);
t57 = t25 * t35;
t55 = t33 * t24;
t54 = t33 * t25;
t37 = -pkin(9) - pkin(8);
t51 = t35 * t37;
t50 = t36 * t24;
t49 = t36 * t25;
t46 = t36 * pkin(1) + t33 * qJ(2);
t44 = pkin(4) * t31 + pkin(7);
t11 = t32 * t50 + t54;
t9 = t32 * t55 - t49;
t43 = g(1) * t11 + g(2) * t9;
t20 = g(1) * t36 + g(2) * t33;
t23 = t34 * pkin(4) + pkin(3);
t42 = t32 * t23 + t51;
t41 = pkin(5) * t25 + qJ(6) * t24 + t23;
t1 = g(1) * t9 - g(2) * t11 + g(3) * t58;
t10 = t32 * t54 + t50;
t12 = t32 * t49 - t55;
t3 = g(1) * t10 - g(2) * t12 + g(3) * t57;
t38 = -g(1) * (-t9 * pkin(5) + t10 * qJ(6)) - g(2) * (t11 * pkin(5) - t12 * qJ(6)) - g(3) * (-pkin(5) * t58 + qJ(6) * t57);
t27 = t36 * qJ(2);
t18 = t20 * t35;
t17 = t32 * t47 - t53;
t15 = t32 * t52 + t48;
t13 = -t32 * t67 + t59;
t6 = t40 * t25;
t5 = t40 * t24;
t4 = -g(1) * t12 - g(2) * t10;
t2 = [0, -t67, t20, t67, -t20, -g(1) * (-t33 * pkin(1) + t27) - g(2) * t46, 0, 0, 0, 0, 0, -t20 * t32, -t18, 0, 0, 0, 0, 0, -g(1) * t17 - g(2) * t15, g(1) * t16 - g(2) * t14, 0, 0, 0, 0, 0, t4, t43, t4, t18, -t43, -g(1) * (t12 * pkin(5) + t11 * qJ(6) + t27) - g(2) * (t10 * pkin(5) + t9 * qJ(6) + t46) + (-g(1) * t42 - g(2) * t44) * t36 + (-g(1) * (-pkin(1) - t44) - g(2) * t42) * t33; 0, 0, 0, 0, 0, t67, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, t13, 0, 0, 0, 0, 0, -t40 * t34, t40 * t31, 0, 0, 0, 0, 0, -t6, t5, -t6, -t13, -t5, -g(3) * (-t41 * t32 - t51) + t67 * (-t32 * t37 + t41 * t35); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, g(1) * t15 - g(2) * t17 + t34 * t59, 0, 0, 0, 0, 0, t1, t3, t1, 0, -t3, t66 * pkin(4) + t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t3, t1, 0, -t3, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t2;
