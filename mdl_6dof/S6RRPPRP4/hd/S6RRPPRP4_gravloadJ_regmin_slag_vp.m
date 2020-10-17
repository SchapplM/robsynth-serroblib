% Calculate minimal parameter regressor of gravitation load for
% S6RRPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
% 
% Output:
% taug_reg [6x29]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPPRP4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:21:49
% EndTime: 2019-05-06 09:21:51
% DurationCPUTime: 0.39s
% Computational Cost: add. (267->99), mult. (704->137), div. (0->0), fcn. (811->8), ass. (0->68)
t54 = sin(qJ(2));
t51 = sin(pkin(9));
t52 = cos(pkin(9));
t53 = sin(qJ(5));
t56 = cos(qJ(5));
t67 = t51 * t53 + t52 * t56;
t95 = t67 * t54;
t55 = sin(qJ(1));
t58 = cos(qJ(1));
t69 = g(1) * t58 + g(2) * t55;
t94 = t69 * t54;
t88 = t51 * t56;
t93 = t54 * (t52 * t53 - t88);
t92 = g(1) * t55;
t89 = g(3) * t54;
t57 = cos(qJ(2));
t47 = t57 * pkin(2);
t86 = t52 * t57;
t85 = t54 * t55;
t84 = t54 * t58;
t83 = t55 * t57;
t82 = t57 * t58;
t81 = t58 * t51;
t80 = t58 * t52;
t45 = t54 * qJ(3);
t79 = t47 + t45;
t78 = qJ(3) * t57;
t77 = qJ(4) * t51;
t76 = -pkin(1) - t47;
t75 = -pkin(2) - t77;
t74 = pkin(3) * t86 + t57 * t77 + t79;
t73 = pkin(2) * t82 + t55 * pkin(7) + (pkin(1) + t45) * t58;
t29 = t51 * t83 + t80;
t30 = t52 * t83 - t81;
t7 = t29 * t56 - t30 * t53;
t31 = -t55 * t52 + t57 * t81;
t32 = t55 * t51 + t57 * t80;
t9 = -t31 * t56 + t32 * t53;
t72 = g(1) * t7 + g(2) * t9;
t48 = t58 * pkin(7);
t71 = -t30 * pkin(3) - t29 * qJ(4) + t48;
t70 = g(1) * t29 - g(2) * t31;
t68 = -g(2) * t58 + t92;
t6 = t29 * t53 + t30 * t56;
t1 = g(1) * t9 - g(2) * t7 + g(3) * t93;
t10 = t31 * t53 + t32 * t56;
t63 = g(1) * t10 + g(2) * t6 + g(3) * t95;
t14 = t55 * t93;
t16 = t58 * t93;
t25 = t53 * t86 - t57 * t88;
t62 = g(1) * t16 + g(2) * t14 - g(3) * t25;
t60 = t32 * pkin(3) + t31 * qJ(4) + t73;
t59 = (t76 - t45) * t92;
t21 = -g(3) * t57 + t94;
t40 = t58 * t78;
t37 = t55 * t78;
t33 = g(1) * t85 - g(2) * t84;
t26 = t67 * t57;
t22 = t69 * t57 + t89;
t17 = t58 * t95;
t15 = t55 * t95;
t13 = t21 * t52;
t12 = t21 * t51;
t11 = g(1) * t30 - g(2) * t32;
t4 = -g(1) * t31 - g(2) * t29 - t51 * t89;
t3 = g(1) * t17 + g(2) * t15 - g(3) * t26;
t2 = g(1) * t6 - g(2) * t10;
t5 = [0, t68, t69, 0, 0, 0, 0, 0, t68 * t57, -t33, t11, -t70, t33, -g(1) * t48 - g(2) * t73 - t59, t11, t33, t70, -g(1) * t71 - g(2) * t60 - t59, 0, 0, 0, 0, 0, t2, t72, t2, -t33, -t72, -g(1) * (-t30 * pkin(4) - pkin(5) * t6 + t7 * qJ(6) + t71) - g(2) * (t32 * pkin(4) + t10 * pkin(5) - pkin(8) * t84 + t9 * qJ(6) + t60) - ((pkin(8) - qJ(3)) * t54 + t76) * t92; 0, 0, 0, 0, 0, 0, 0, 0, t21, t22, t13, -t12, -t22, -g(1) * (-pkin(2) * t84 + t40) - g(2) * (-pkin(2) * t85 + t37) - g(3) * t79, t13, -t22, t12, -g(1) * t40 - g(2) * t37 - g(3) * t74 + (pkin(3) * t52 - t75) * t94, 0, 0, 0, 0, 0, t3, -t62, t3, t22, t62, -g(1) * (-t17 * pkin(5) - pkin(8) * t82 - t16 * qJ(6) + t40) - g(2) * (-t15 * pkin(5) - pkin(8) * t83 - t14 * qJ(6) + t37) - g(3) * (pkin(4) * t86 + t26 * pkin(5) + t25 * qJ(6) + t74) + (g(3) * pkin(8) + t69 * (-(-pkin(3) - pkin(4)) * t52 - t75)) * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, 0, 0, 0, -t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t63, t1, 0, -t63, -g(1) * (-t9 * pkin(5) + t10 * qJ(6)) - g(2) * (pkin(5) * t7 + t6 * qJ(6)) - g(3) * (-pkin(5) * t93 + qJ(6) * t95); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t5;
