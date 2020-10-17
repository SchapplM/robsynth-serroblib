% Calculate minimal parameter regressor of gravitation load for
% S6RRRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
% 
% Output:
% taug_reg [6x35]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 20:30:53
% EndTime: 2019-05-07 20:30:54
% DurationCPUTime: 0.44s
% Computational Cost: add. (359->70), mult. (532->102), div. (0->0), fcn. (599->10), ass. (0->55)
t32 = qJ(2) + qJ(3);
t30 = cos(t32);
t34 = sin(qJ(4));
t40 = cos(qJ(1));
t57 = t40 * t34;
t36 = sin(qJ(1));
t38 = cos(qJ(4));
t58 = t36 * t38;
t17 = t30 * t57 - t58;
t56 = t40 * t38;
t59 = t36 * t34;
t18 = t30 * t56 + t59;
t33 = sin(qJ(6));
t37 = cos(qJ(6));
t4 = t17 * t37 - t18 * t33;
t48 = t33 * t38 - t34 * t37;
t29 = sin(t32);
t63 = g(3) * t29;
t15 = t30 * t59 + t56;
t16 = t30 * t58 - t57;
t76 = t15 * t37 - t16 * t33;
t80 = g(1) * t4 + g(2) * t76 - t48 * t63;
t79 = pkin(4) * t38 + qJ(5) * t34;
t51 = g(1) * t40 + g(2) * t36;
t78 = t51 * t29;
t9 = -g(3) * t30 + t78;
t77 = t30 * pkin(3) + t29 * pkin(9);
t47 = t33 * t34 + t37 * t38;
t49 = t15 * t33 + t16 * t37;
t5 = t17 * t33 + t18 * t37;
t73 = g(1) * t5 + g(2) * t49 + t47 * t63;
t35 = sin(qJ(2));
t70 = pkin(2) * t35;
t68 = pkin(9) * t30;
t61 = g(3) * t34;
t53 = t79 * t30 + t77;
t52 = g(1) * t15 - g(2) * t17;
t50 = g(1) * t36 - g(2) * t40;
t39 = cos(qJ(2));
t31 = t39 * pkin(2);
t46 = t31 + pkin(1) + t77;
t3 = g(1) * t17 + g(2) * t15 + t29 * t61;
t43 = g(1) * t18 + g(2) * t16 + t38 * t63;
t42 = (pkin(3) + t79) * t78;
t41 = -pkin(8) - pkin(7);
t22 = t40 * t68;
t20 = t36 * t68;
t14 = t50 * t29;
t10 = t51 * t30 + t63;
t8 = t9 * t38;
t7 = -t30 * t61 + t34 * t78;
t6 = g(1) * t16 - g(2) * t18;
t2 = t9 * t47;
t1 = t9 * t48;
t11 = [0, t50, t51, 0, 0, 0, 0, 0, t50 * t39, -t50 * t35, 0, 0, 0, 0, 0, t50 * t30, -t14, 0, 0, 0, 0, 0, t6, -t52, t6, t14, t52, -g(1) * (-t16 * pkin(4) - t15 * qJ(5)) - g(2) * (t18 * pkin(4) + t17 * qJ(5)) + (g(1) * t41 - g(2) * t46) * t40 + (g(1) * t46 + g(2) * t41) * t36, 0, 0, 0, 0, 0, g(1) * t49 - g(2) * t5, g(1) * t76 - g(2) * t4; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t39 + t51 * t35, g(3) * t35 + t51 * t39, 0, 0, 0, 0, 0, t9, t10, 0, 0, 0, 0, 0, t8, -t7, t8, -t10, t7, -g(1) * (-t40 * t70 + t22) - g(2) * (-t36 * t70 + t20) - g(3) * (t31 + t53) + t42, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, 0, 0, 0, 0, t8, -t7, t8, -t10, t7, -g(1) * t22 - g(2) * t20 - g(3) * t53 + t42, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t43, t3, 0, -t43, -g(1) * (-t17 * pkin(4) + t18 * qJ(5)) - g(2) * (-t15 * pkin(4) + t16 * qJ(5)) - (-pkin(4) * t34 + qJ(5) * t38) * t63, 0, 0, 0, 0, 0, t80, -t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80, t73;];
taug_reg  = t11;
