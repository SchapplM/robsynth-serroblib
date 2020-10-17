% Calculate minimal parameter regressor of gravitation load for
% S6PRRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% 
% Output:
% taug_reg [6x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRPRP5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 04:10:55
% EndTime: 2019-05-05 04:10:56
% DurationCPUTime: 0.37s
% Computational Cost: add. (355->98), mult. (940->143), div. (0->0), fcn. (1163->10), ass. (0->66)
t57 = sin(qJ(3));
t91 = qJ(4) * t57 + pkin(2);
t54 = sin(pkin(10));
t58 = sin(qJ(2));
t61 = cos(qJ(2));
t75 = cos(pkin(10));
t76 = cos(pkin(6));
t68 = t76 * t75;
t39 = t54 * t61 + t58 * t68;
t60 = cos(qJ(3));
t55 = sin(pkin(6));
t70 = t55 * t75;
t22 = t39 * t60 - t57 * t70;
t71 = t54 * t76;
t41 = -t58 * t71 + t75 * t61;
t24 = t54 * t55 * t57 + t41 * t60;
t83 = t55 * t60;
t43 = t76 * t57 + t58 * t83;
t64 = g(1) * t24 + g(2) * t22 + g(3) * t43;
t90 = pkin(4) + pkin(8);
t38 = t54 * t58 - t61 * t68;
t86 = t38 * t60;
t40 = t75 * t58 + t61 * t71;
t85 = t40 * t60;
t84 = t55 * t58;
t82 = t55 * t61;
t56 = sin(qJ(5));
t81 = t56 * t57;
t80 = t56 * t61;
t59 = cos(qJ(5));
t79 = t57 * t59;
t78 = t60 * t61;
t74 = t59 * t82;
t73 = -pkin(3) * t86 - t38 * t91;
t72 = -pkin(3) * t85 - t40 * t91;
t69 = t55 * pkin(3) * t78 + pkin(8) * t84 + t82 * t91;
t42 = t57 * t84 - t76 * t60;
t25 = t42 * t59 + t55 * t80;
t21 = t39 * t57 + t60 * t70;
t7 = -t21 * t59 + t38 * t56;
t23 = t41 * t57 - t54 * t83;
t9 = -t23 * t59 + t40 * t56;
t1 = g(1) * t9 + g(2) * t7 - g(3) * t25;
t10 = t23 * t56 + t40 * t59;
t26 = -t42 * t56 + t74;
t8 = t21 * t56 + t38 * t59;
t66 = g(1) * t10 + g(2) * t8 - g(3) * t26;
t13 = t38 * t79 + t39 * t56;
t15 = t40 * t79 + t41 * t56;
t29 = t56 * t84 - t57 * t74;
t65 = g(1) * t15 + g(2) * t13 + g(3) * t29;
t6 = g(1) * t23 + g(2) * t21 + g(3) * t42;
t63 = -g(1) * t40 - g(2) * t38 + g(3) * t82;
t62 = g(1) * t41 + g(2) * t39 + g(3) * t84;
t37 = t42 * pkin(3);
t30 = (t57 * t80 + t58 * t59) * t55;
t20 = t23 * pkin(3);
t19 = t21 * pkin(3);
t16 = -t40 * t81 + t41 * t59;
t14 = -t38 * t81 + t39 * t59;
t12 = t63 * t60;
t11 = t63 * t57;
t4 = t64 * t59;
t3 = t64 * t56;
t2 = -g(1) * t16 - g(2) * t14 - g(3) * t30;
t5 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, -t63, t62, 0, 0, 0, 0, 0, -t12, t11, -t62, t12, -t11, -g(1) * (t41 * pkin(8) + t72) - g(2) * (t39 * pkin(8) + t73) - g(3) * t69, 0, 0, 0, 0, 0, t2, t65, t2, -t12, -t65, -g(1) * (t16 * pkin(5) - pkin(9) * t85 + t15 * qJ(6) + t90 * t41 + t72) - g(2) * (t14 * pkin(5) - pkin(9) * t86 + t13 * qJ(6) + t90 * t39 + t73) - g(3) * (t30 * pkin(5) + t29 * qJ(6) + (pkin(4) * t58 + pkin(9) * t78) * t55 + t69); 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t64, 0, -t6, -t64, -g(1) * (t24 * qJ(4) - t20) - g(2) * (t22 * qJ(4) - t19) - g(3) * (t43 * qJ(4) - t37) 0, 0, 0, 0, 0, -t3, -t4, -t3, t6, t4, -g(1) * (-t23 * pkin(9) - t20) - g(2) * (-t21 * pkin(9) - t19) - g(3) * (-t42 * pkin(9) - t37) - t64 * (pkin(5) * t56 - qJ(6) * t59 + qJ(4)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t66, t1, 0, -t66, -g(1) * (-t9 * pkin(5) + t10 * qJ(6)) - g(2) * (-t7 * pkin(5) + t8 * qJ(6)) - g(3) * (t25 * pkin(5) - t26 * qJ(6)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t5;
