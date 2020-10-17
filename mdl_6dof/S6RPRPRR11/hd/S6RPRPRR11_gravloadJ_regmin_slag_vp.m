% Calculate minimal parameter regressor of gravitation load for
% S6RPRPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRR11_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR11_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_gravloadJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 20:18:45
% EndTime: 2019-05-05 20:18:47
% DurationCPUTime: 0.72s
% Computational Cost: add. (567->94), mult. (1425->170), div. (0->0), fcn. (1850->16), ass. (0->60)
t45 = cos(qJ(1));
t69 = sin(pkin(12));
t73 = cos(pkin(6));
t60 = t73 * t69;
t71 = cos(pkin(12));
t78 = sin(qJ(1));
t26 = t45 * t60 + t71 * t78;
t43 = sin(qJ(3));
t79 = cos(qJ(3));
t61 = t73 * t71;
t51 = -t45 * t61 + t78 * t69;
t40 = sin(pkin(6));
t70 = sin(pkin(7));
t67 = t40 * t70;
t72 = cos(pkin(7));
t85 = t45 * t67 + t51 * t72;
t10 = t26 * t43 + t85 * t79;
t42 = sin(qJ(6));
t44 = cos(qJ(6));
t13 = -t26 * t79 + t85 * t43;
t68 = t40 * t72;
t19 = -t45 * t68 + t51 * t70;
t38 = pkin(13) + qJ(5);
t35 = sin(t38);
t36 = cos(t38);
t5 = t13 * t36 - t19 * t35;
t89 = t10 * t44 + t5 * t42;
t88 = -t10 * t42 + t5 * t44;
t84 = t13 * t35 + t19 * t36;
t47 = t45 * t69 + t61 * t78;
t80 = t47 * t72 - t78 * t67;
t77 = t36 * t42;
t76 = t36 * t44;
t74 = qJ(2) * t40;
t75 = t45 * pkin(1) + t78 * t74;
t64 = -pkin(1) * t78 + t45 * t74;
t27 = t45 * t71 - t60 * t78;
t14 = t27 * t43 + t80 * t79;
t63 = -g(1) * t10 + g(2) * t14;
t59 = t72 * t71;
t58 = t70 * t73;
t56 = g(1) * t78 - g(2) * t45;
t55 = -g(1) * t45 - g(2) * t78;
t17 = t43 * t58 + (t43 * t59 + t69 * t79) * t40;
t25 = -t67 * t71 + t72 * t73;
t15 = t27 * t79 - t80 * t43;
t20 = t47 * t70 + t68 * t78;
t6 = -t15 * t35 + t20 * t36;
t54 = g(1) * t6 + g(2) * t84 + g(3) * (-t17 * t35 + t25 * t36);
t16 = -t79 * t58 + (t43 * t69 - t59 * t79) * t40;
t53 = g(1) * t14 + g(2) * t10 + g(3) * t16;
t52 = g(1) * t15 - g(2) * t13 + g(3) * t17;
t41 = cos(pkin(13));
t39 = sin(pkin(13));
t24 = -g(3) * t73 - t40 * t56;
t9 = t17 * t36 + t25 * t35;
t7 = t15 * t36 + t20 * t35;
t2 = t14 * t42 + t7 * t44;
t1 = t14 * t44 - t7 * t42;
t3 = [0, t56, -t55, g(1) * t26 - g(2) * t27, -g(1) * t51 + g(2) * t47, t55 * t40, -g(1) * t64 - g(2) * t75, 0, 0, 0, 0, 0, -g(1) * t13 - g(2) * t15, t63, -g(1) * (t13 * t41 - t19 * t39) - g(2) * (t15 * t41 + t20 * t39) -g(1) * (-t13 * t39 - t19 * t41) - g(2) * (-t15 * t39 + t20 * t41) -t63, -g(1) * (-t26 * pkin(2) + t13 * pkin(3) - qJ(4) * t10 + t64) - g(2) * (t27 * pkin(2) + t15 * pkin(3) + t14 * qJ(4) + t75) + (g(1) * t19 - g(2) * t20) * pkin(9), 0, 0, 0, 0, 0, -g(1) * t5 - g(2) * t7, g(1) * t84 - g(2) * t6, 0, 0, 0, 0, 0, -g(1) * t88 - g(2) * t2, g(1) * t89 - g(2) * t1; 0, 0, 0, 0, 0, 0, t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, t52, t53 * t41, -t53 * t39, -t52, -g(1) * (-t14 * pkin(3) + t15 * qJ(4)) - g(2) * (-t10 * pkin(3) - qJ(4) * t13) - g(3) * (-t16 * pkin(3) + t17 * qJ(4)) 0, 0, 0, 0, 0, t53 * t36, -t53 * t35, 0, 0, 0, 0, 0, -g(1) * (-t14 * t76 + t15 * t42) - g(2) * (-t10 * t76 - t13 * t42) - g(3) * (-t16 * t76 + t17 * t42) -g(1) * (t14 * t77 + t15 * t44) - g(2) * (t10 * t77 - t13 * t44) - g(3) * (t16 * t77 + t17 * t44); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, g(1) * t7 - g(2) * t5 + g(3) * t9, 0, 0, 0, 0, 0, -t54 * t44, t54 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 - g(2) * t89 - g(3) * (t16 * t44 - t9 * t42) g(1) * t2 - g(2) * t88 - g(3) * (-t16 * t42 - t9 * t44);];
taug_reg  = t3;
