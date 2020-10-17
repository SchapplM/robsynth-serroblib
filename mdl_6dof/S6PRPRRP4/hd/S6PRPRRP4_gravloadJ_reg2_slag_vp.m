% Calculate inertial parameters regressor of gravitation load for
% S6PRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRRP4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:50:28
% EndTime: 2019-05-04 23:50:30
% DurationCPUTime: 0.44s
% Computational Cost: add. (547->103), mult. (963->147), div. (0->0), fcn. (1176->12), ass. (0->69)
t57 = pkin(11) + qJ(4);
t55 = sin(t57);
t56 = cos(t57);
t96 = pkin(4) * t56 + pkin(9) * t55;
t63 = sin(qJ(5));
t93 = t56 * t63;
t65 = cos(qJ(5));
t92 = t56 * t65;
t59 = sin(pkin(10));
t60 = sin(pkin(6));
t91 = t59 * t60;
t64 = sin(qJ(2));
t90 = t60 * t64;
t66 = cos(qJ(2));
t89 = t60 * t66;
t88 = t65 * t66;
t84 = cos(pkin(10));
t85 = cos(pkin(6));
t73 = t85 * t84;
t42 = t59 * t64 - t66 * t73;
t43 = t59 * t66 + t64 * t73;
t61 = cos(pkin(11));
t54 = t61 * pkin(3) + pkin(2);
t62 = -pkin(8) - qJ(3);
t87 = -t42 * t54 - t43 * t62;
t79 = t59 * t85;
t44 = t84 * t64 + t66 * t79;
t45 = -t64 * t79 + t84 * t66;
t86 = -t44 * t54 - t45 * t62;
t83 = t63 * t89;
t78 = t60 * t84;
t19 = -t43 * t55 - t56 * t78;
t20 = t43 * t56 - t55 * t78;
t82 = t19 * pkin(4) + t20 * pkin(9);
t21 = -t45 * t55 + t56 * t91;
t22 = t45 * t56 + t55 * t91;
t81 = t21 * pkin(4) + t22 * pkin(9);
t34 = -t55 * t90 + t85 * t56;
t35 = t85 * t55 + t56 * t90;
t80 = t34 * pkin(4) + t35 * pkin(9);
t77 = -t96 * t42 + t87;
t76 = -t96 * t44 + t86;
t75 = t54 * t89 - t62 * t90;
t74 = pkin(5) * t65 + qJ(6) * t63;
t72 = t96 * t89 + t75;
t23 = t35 * t63 + t60 * t88;
t7 = t20 * t63 - t42 * t65;
t9 = t22 * t63 - t44 * t65;
t1 = g(1) * t9 + g(2) * t7 + g(3) * t23;
t10 = t22 * t65 + t44 * t63;
t24 = t35 * t65 - t83;
t8 = t20 * t65 + t42 * t63;
t71 = g(1) * t10 + g(2) * t8 + g(3) * t24;
t11 = -t42 * t93 - t43 * t65;
t13 = -t44 * t93 - t45 * t65;
t25 = t56 * t83 - t65 * t90;
t70 = g(1) * t13 + g(2) * t11 + g(3) * t25;
t69 = g(1) * t21 + g(2) * t19 + g(3) * t34;
t68 = g(1) * t22 + g(2) * t20 + g(3) * t35;
t15 = -g(1) * t44 - g(2) * t42 + g(3) * t89;
t67 = g(1) * t45 + g(2) * t43 + g(3) * t90;
t26 = (t56 * t88 + t63 * t64) * t60;
t14 = -t44 * t92 + t45 * t63;
t12 = -t42 * t92 + t43 * t63;
t6 = t15 * t55;
t4 = t69 * t65;
t3 = t69 * t63;
t2 = -g(1) * t14 - g(2) * t12 - g(3) * t26;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, t67, 0, 0, 0, 0, 0, 0, 0, 0, -t15 * t61, t15 * sin(pkin(11)) -t67, -g(1) * (-t44 * pkin(2) + t45 * qJ(3)) - g(2) * (-t42 * pkin(2) + t43 * qJ(3)) - g(3) * (pkin(2) * t66 + qJ(3) * t64) * t60, 0, 0, 0, 0, 0, 0, -t15 * t56, t6, -t67, -g(1) * t86 - g(2) * t87 - g(3) * t75, 0, 0, 0, 0, 0, 0, t2, t70, -t6, -g(1) * t76 - g(2) * t77 - g(3) * t72, 0, 0, 0, 0, 0, 0, t2, -t6, -t70, -g(1) * (t14 * pkin(5) + t13 * qJ(6) + t76) - g(2) * (t12 * pkin(5) + t11 * qJ(6) + t77) - g(3) * (t26 * pkin(5) + t25 * qJ(6) + t72); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, t68, 0, 0, 0, 0, 0, 0, 0, 0, -t4, t3, -t68, -g(1) * t81 - g(2) * t82 - g(3) * t80, 0, 0, 0, 0, 0, 0, -t4, -t68, -t3, -g(1) * (t74 * t21 + t81) - g(2) * (t74 * t19 + t82) - g(3) * (t74 * t34 + t80); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t71, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t71, -g(1) * (-t9 * pkin(5) + t10 * qJ(6)) - g(2) * (-t7 * pkin(5) + t8 * qJ(6)) - g(3) * (-t23 * pkin(5) + t24 * qJ(6)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t5;
