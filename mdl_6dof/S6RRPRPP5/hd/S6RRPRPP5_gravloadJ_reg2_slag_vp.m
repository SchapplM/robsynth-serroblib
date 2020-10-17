% Calculate inertial parameters regressor of gravitation load for
% S6RRPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRPP5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:44:59
% EndTime: 2019-05-06 12:45:00
% DurationCPUTime: 0.40s
% Computational Cost: add. (232->101), mult. (564->119), div. (0->0), fcn. (584->6), ass. (0->61)
t42 = sin(qJ(1));
t45 = cos(qJ(1));
t22 = g(1) * t45 + g(2) * t42;
t41 = sin(qJ(2));
t44 = cos(qJ(2));
t9 = g(3) * t41 + t22 * t44;
t83 = pkin(2) + pkin(8);
t82 = pkin(2) * t41;
t81 = g(1) * t42;
t77 = g(3) * t44;
t35 = t44 * pkin(2);
t40 = sin(qJ(4));
t76 = t40 * t44;
t75 = t42 * t40;
t43 = cos(qJ(4));
t74 = t42 * t43;
t73 = t44 * t45;
t72 = t45 * t40;
t71 = t45 * t43;
t65 = qJ(3) * t44;
t24 = t42 * t65;
t62 = pkin(4) * t76;
t70 = t42 * t62 + t24;
t26 = t45 * t65;
t69 = t45 * t62 + t26;
t31 = t41 * qJ(3);
t68 = t35 + t31;
t36 = t45 * pkin(7);
t67 = t45 * pkin(3) + t36;
t66 = t45 * pkin(1) + t42 * pkin(7);
t64 = qJ(5) * t40;
t63 = qJ(5) * t43;
t61 = -qJ(6) + t83;
t60 = t44 * pkin(8) + t68;
t59 = t44 * t63;
t58 = -pkin(1) - t31;
t14 = -t41 * t71 + t75;
t15 = t41 * t72 + t74;
t57 = -t14 * pkin(4) + t15 * qJ(5);
t16 = t41 * t74 + t72;
t17 = -t41 * t75 + t71;
t56 = t16 * pkin(4) - t17 * qJ(5);
t55 = pkin(2) * t73 + t45 * t31 + t66;
t54 = g(1) * t16 + g(2) * t14;
t53 = -g(2) * t45 + t81;
t52 = pkin(5) * t40 - t63;
t51 = g(3) * (t41 * t40 * pkin(4) + t60);
t50 = t17 * pkin(4) + t16 * qJ(5) + t67;
t49 = t42 * pkin(3) + pkin(8) * t73 + t55;
t2 = g(1) * t14 - g(2) * t16 + t43 * t77;
t3 = -g(1) * t15 + g(2) * t17 + g(3) * t76;
t48 = t22 * t83;
t47 = t15 * pkin(4) + t14 * qJ(5) + t49;
t46 = (-t83 * t44 + t58) * t81;
t19 = -g(2) * t73 + t44 * t81;
t18 = t53 * t41;
t8 = t22 * t41 - t77;
t7 = t9 * t43;
t6 = t9 * t40;
t5 = -g(1) * t17 - g(2) * t15;
t1 = [0, 0, 0, 0, 0, 0, t53, t22, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t18, -t22, -g(1) * (-t42 * pkin(1) + t36) - g(2) * t66, 0, 0, 0, 0, 0, 0, -t22, -t19, t18, -g(1) * t36 - g(2) * t55 - (t58 - t35) * t81, 0, 0, 0, 0, 0, 0, t5, t54, t19, -g(1) * t67 - g(2) * t49 - t46, 0, 0, 0, 0, 0, 0, t5, t19, -t54, -g(1) * t50 - g(2) * t47 - t46, 0, 0, 0, 0, 0, 0, t5, -t54, -t19, -g(1) * (t17 * pkin(5) + t50) - g(2) * (t15 * pkin(5) - qJ(6) * t73 + t47) - (-t61 * t44 + t58) * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, -t9, -g(1) * (-t45 * t82 + t26) - g(2) * (-t42 * t82 + t24) - g(3) * t68, 0, 0, 0, 0, 0, 0, -t6, -t7, t8, -g(1) * t26 - g(2) * t24 - g(3) * t60 + t48 * t41, 0, 0, 0, 0, 0, 0, -t6, t8, t7, -g(1) * (-t45 * t59 + t69) - g(2) * (-t42 * t59 + t70) - t51 + (g(3) * t63 + t48) * t41, 0, 0, 0, 0, 0, 0, -t6, t7, -t8, -g(1) * t69 - g(2) * t70 - t51 + (g(3) * qJ(6) - t22 * t52) * t44 + (-g(3) * t52 + t22 * t61) * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t3, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, t3, -g(1) * t57 - g(2) * t56 - (-pkin(4) * t43 - t64) * t77, 0, 0, 0, 0, 0, 0, t2, t3, 0, -g(1) * (-t14 * pkin(5) + t57) - g(2) * (t16 * pkin(5) + t56) - (-t64 + (-pkin(4) - pkin(5)) * t43) * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9;];
taug_reg  = t1;
