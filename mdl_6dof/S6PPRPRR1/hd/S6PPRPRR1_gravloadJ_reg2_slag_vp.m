% Calculate inertial parameters regressor of gravitation load for
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
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PPRPRR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRPRR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_gravloadJ_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:53:55
% EndTime: 2019-05-04 19:53:57
% DurationCPUTime: 0.56s
% Computational Cost: add. (747->108), mult. (2073->184), div. (0->0), fcn. (2705->16), ass. (0->73)
t46 = sin(pkin(7));
t43 = sin(pkin(13));
t48 = cos(pkin(13));
t55 = sin(qJ(3));
t58 = cos(qJ(3));
t73 = t43 * t58 + t48 * t55;
t29 = t73 * t46;
t51 = cos(pkin(7));
t31 = t73 * t51;
t37 = t43 * t55 - t48 * t58;
t44 = sin(pkin(12));
t47 = sin(pkin(6));
t49 = cos(pkin(12));
t52 = cos(pkin(6));
t20 = t52 * t29 + (t31 * t49 - t37 * t44) * t47;
t45 = sin(pkin(11));
t50 = cos(pkin(11));
t82 = t50 * t52;
t33 = -t44 * t45 + t49 * t82;
t87 = t45 * t52;
t35 = -t44 * t50 - t49 * t87;
t83 = t47 * t51;
t76 = t49 * t83;
t84 = t47 * t50;
t86 = t46 * t47;
t92 = g(1) * (t35 * t51 + t45 * t86) - g(2) * (-t33 * t51 + t46 * t84) + g(3) * (t46 * t52 + t76);
t91 = pkin(3) * t55;
t90 = pkin(3) * t58;
t89 = t44 * t47;
t88 = t45 * t47;
t85 = t46 * t58;
t53 = sin(qJ(6));
t57 = cos(qJ(5));
t80 = t53 * t57;
t56 = cos(qJ(6));
t79 = t56 * t57;
t78 = t51 * t90;
t77 = t47 * t85;
t54 = sin(qJ(5));
t75 = pkin(5) * t57 + pkin(10) * t54;
t36 = -t44 * t87 + t49 * t50;
t72 = pkin(3) * t45 * t77 + t35 * t78 - t36 * t91;
t68 = pkin(3) * t52 * t85 + t76 * t90 - t89 * t91;
t32 = -t49 * t86 + t51 * t52;
t14 = -t20 * t54 + t32 * t57;
t21 = -t33 * t46 - t50 * t83;
t34 = t44 * t82 + t45 * t49;
t8 = t29 * t84 - t31 * t33 + t34 * t37;
t2 = t21 * t57 + t54 * t8;
t13 = t29 * t88 + t31 * t35 - t36 * t37;
t22 = -t35 * t46 + t45 * t83;
t4 = -t13 * t54 + t22 * t57;
t67 = g(1) * t4 + g(2) * t2 + g(3) * t14;
t15 = t20 * t57 + t32 * t54;
t3 = t21 * t54 - t57 * t8;
t5 = t13 * t57 + t22 * t54;
t66 = g(1) * t5 + g(2) * t3 + g(3) * t15;
t65 = -g(1) * t13 + g(2) * t8 - g(3) * t20;
t28 = t37 * t46;
t30 = t37 * t51;
t12 = -t28 * t88 - t30 * t35 - t36 * t73;
t19 = -t52 * t28 + (-t30 * t49 - t44 * t73) * t47;
t9 = t28 * t84 - t30 * t33 - t34 * t73;
t64 = g(1) * t12 + g(2) * t9 + g(3) * t19;
t63 = g(1) * t36 + g(2) * t34 + g(3) * t89;
t62 = pkin(4) * t12 + pkin(9) * t13 + t72;
t61 = pkin(4) * t19 + pkin(9) * t20 + t68;
t60 = t33 * t78 + (-t34 * t55 - t50 * t77) * pkin(3);
t59 = pkin(4) * t9 - t8 * pkin(9) + t60;
t27 = -g(3) * t52 + (-g(1) * t45 + g(2) * t50) * t47;
t16 = -g(1) * t22 - g(2) * t21 - g(3) * t32;
t1 = t64 * t54;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55 * t63 - t58 * t92, t55 * t92 + t58 * t63, 0, 0, 0, 0, 0, 0, 0, 0, -t64, -t65, 0, -g(1) * t72 - g(2) * t60 - g(3) * t68, 0, 0, 0, 0, 0, 0, -t64 * t57, t1, t65, -g(1) * t62 - g(2) * t59 - g(3) * t61, 0, 0, 0, 0, 0, 0, -g(1) * (t12 * t79 + t13 * t53) - g(2) * (-t53 * t8 + t79 * t9) - g(3) * (t19 * t79 + t20 * t53) -g(1) * (-t12 * t80 + t13 * t56) - g(2) * (-t56 * t8 - t80 * t9) - g(3) * (-t19 * t80 + t20 * t56) -t1, -g(1) * (t12 * t75 + t62) - g(2) * (t75 * t9 + t59) - g(3) * (t19 * t75 + t61); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, t66, 0, 0, 0, 0, 0, 0, 0, 0, -t67 * t56, t67 * t53, -t66, -g(1) * (pkin(5) * t4 + pkin(10) * t5) - g(2) * (pkin(5) * t2 + pkin(10) * t3) - g(3) * (pkin(5) * t14 + pkin(10) * t15); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t12 * t56 - t5 * t53) - g(2) * (-t3 * t53 - t56 * t9) - g(3) * (-t15 * t53 - t19 * t56) -g(1) * (t12 * t53 - t5 * t56) - g(2) * (-t3 * t56 + t53 * t9) - g(3) * (-t15 * t56 + t19 * t53) 0, 0;];
taug_reg  = t6;
