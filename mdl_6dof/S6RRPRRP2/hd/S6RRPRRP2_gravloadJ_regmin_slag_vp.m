% Calculate minimal parameter regressor of gravitation load for
% S6RRPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% taug_reg [6x30]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRP2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP2_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP2_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 17:25:07
% EndTime: 2019-05-06 17:25:08
% DurationCPUTime: 0.25s
% Computational Cost: add. (366->65), mult. (355->87), div. (0->0), fcn. (365->10), ass. (0->48)
t35 = sin(qJ(5));
t38 = cos(qJ(5));
t62 = pkin(5) * t38 + qJ(6) * t35;
t37 = sin(qJ(1));
t40 = cos(qJ(1));
t20 = g(1) * t40 + g(2) * t37;
t33 = qJ(2) + pkin(10);
t29 = qJ(4) + t33;
t24 = sin(t29);
t44 = t20 * t24;
t25 = cos(t29);
t61 = t25 * pkin(4) + t24 * pkin(9);
t59 = pkin(9) * t25;
t56 = g(3) * t24;
t55 = g(3) * t35;
t54 = t37 * t35;
t53 = t37 * t38;
t52 = t40 * t35;
t51 = t40 * t38;
t34 = -qJ(3) - pkin(7);
t39 = cos(qJ(2));
t30 = t39 * pkin(2);
t50 = pkin(3) * cos(t33) + t30;
t48 = t62 * t25 + t61;
t10 = t25 * t52 - t53;
t8 = t25 * t54 + t51;
t47 = g(1) * t8 - g(2) * t10;
t19 = g(1) * t37 - g(2) * t40;
t46 = pkin(1) + t50 + t61;
t1 = g(1) * t10 + g(2) * t8 + t24 * t55;
t11 = t25 * t51 + t54;
t9 = t25 * t53 - t52;
t43 = g(1) * t11 + g(2) * t9 + t38 * t56;
t5 = -g(3) * t25 + t44;
t36 = sin(qJ(2));
t42 = -g(3) * t39 + t20 * t36;
t41 = (pkin(4) + t62) * t44;
t32 = -pkin(8) + t34;
t28 = t30 + pkin(1);
t17 = t40 * t59;
t15 = t37 * t59;
t13 = -pkin(3) * sin(t33) - t36 * pkin(2);
t7 = t19 * t24;
t6 = t20 * t25 + t56;
t4 = t5 * t38;
t3 = -t25 * t55 + t35 * t44;
t2 = g(1) * t9 - g(2) * t11;
t12 = [0, t19, t20, 0, 0, 0, 0, 0, t19 * t39, -t19 * t36, -t20, -g(1) * (-t37 * t28 - t40 * t34) - g(2) * (t40 * t28 - t37 * t34) 0, 0, 0, 0, 0, t19 * t25, -t7, 0, 0, 0, 0, 0, t2, -t47, t2, t7, t47, -g(1) * (-t9 * pkin(5) - t8 * qJ(6)) - g(2) * (t11 * pkin(5) + t10 * qJ(6)) + (g(1) * t32 - g(2) * t46) * t40 + (g(1) * t46 + g(2) * t32) * t37; 0, 0, 0, 0, 0, 0, 0, 0, t42, g(3) * t36 + t20 * t39, 0, t42 * pkin(2), 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t4, -t3, t4, -t6, t3, -g(1) * (t40 * t13 + t17) - g(2) * (t37 * t13 + t15) - g(3) * (t48 + t50) + t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t4, -t3, t4, -t6, t3, -g(1) * t17 - g(2) * t15 - g(3) * t48 + t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t43, t1, 0, -t43, -g(1) * (-t10 * pkin(5) + t11 * qJ(6)) - g(2) * (-t8 * pkin(5) + t9 * qJ(6)) - (-pkin(5) * t35 + qJ(6) * t38) * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t12;
