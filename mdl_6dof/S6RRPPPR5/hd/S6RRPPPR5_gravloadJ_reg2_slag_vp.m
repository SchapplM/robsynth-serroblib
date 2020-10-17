% Calculate inertial parameters regressor of gravitation load for
% S6RRPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPPPR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 08:54:16
% EndTime: 2019-05-06 08:54:18
% DurationCPUTime: 0.50s
% Computational Cost: add. (250->103), mult. (631->131), div. (0->0), fcn. (681->8), ass. (0->65)
t45 = cos(qJ(2));
t43 = sin(qJ(1));
t46 = cos(qJ(1));
t18 = g(1) * t46 + g(2) * t43;
t42 = sin(qJ(2));
t87 = t18 * t42;
t9 = -g(3) * t45 + t87;
t86 = g(1) * t43;
t83 = g(3) * t42;
t35 = t45 * pkin(2);
t40 = cos(pkin(9));
t81 = t40 * t45;
t80 = t42 * t43;
t79 = t42 * t46;
t78 = t43 * t45;
t77 = t45 * t46;
t39 = sin(pkin(9));
t76 = t46 * t39;
t75 = t46 * t40;
t74 = -pkin(4) - qJ(3);
t73 = pkin(5) + qJ(4);
t68 = qJ(3) * t45;
t21 = t43 * t68;
t72 = pkin(4) * t78 + t21;
t24 = t46 * t68;
t71 = pkin(4) * t77 + t24;
t32 = t42 * qJ(3);
t70 = t35 + t32;
t69 = t46 * pkin(1) + t43 * pkin(7);
t67 = qJ(4) * t39;
t15 = -t43 * t40 + t45 * t76;
t66 = t15 * qJ(4);
t65 = -pkin(1) - t35;
t64 = pkin(3) * t81 + t45 * t67 + t70;
t63 = pkin(2) * t77 + t46 * t32 + t69;
t13 = t39 * t78 + t75;
t14 = t40 * t78 - t76;
t36 = t46 * pkin(7);
t62 = -t14 * pkin(3) - t13 * qJ(4) + t36;
t16 = t43 * t39 + t45 * t75;
t61 = t16 * pkin(3) + t63;
t5 = g(1) * t13 - g(2) * t15;
t6 = g(1) * t14 - g(2) * t16;
t60 = -g(2) * t46 + t86;
t59 = -pkin(2) + (-pkin(3) - qJ(5)) * t40;
t41 = sin(qJ(6));
t44 = cos(qJ(6));
t58 = t13 * t44 + t14 * t41;
t57 = t13 * t41 - t14 * t44;
t56 = t39 * t44 + t40 * t41;
t55 = t39 * t41 - t40 * t44;
t54 = t42 * pkin(4) + qJ(5) * t81 + t64;
t52 = g(3) * t55;
t51 = -t14 * qJ(5) + t62;
t49 = (t65 - t32) * t86;
t48 = pkin(4) * t79 + t16 * qJ(5) + t61;
t17 = g(1) * t80 - g(2) * t79;
t10 = t18 * t45 + t83;
t8 = t9 * t40;
t7 = t9 * t39;
t4 = t15 * t44 + t16 * t41;
t3 = -t15 * t41 + t16 * t44;
t2 = -g(1) * t16 - g(2) * t14 - t40 * t83;
t1 = -g(1) * t15 - g(2) * t13 - t39 * t83;
t11 = [0, 0, 0, 0, 0, 0, t60, t18, 0, 0, 0, 0, 0, 0, 0, 0, t60 * t45, -t17, -t18, -g(1) * (-t43 * pkin(1) + t36) - g(2) * t69, 0, 0, 0, 0, 0, 0, t6, -t5, t17, -g(1) * t36 - g(2) * t63 - t49, 0, 0, 0, 0, 0, 0, t17, -t6, t5, -g(1) * t62 - g(2) * (t61 + t66) - t49, 0, 0, 0, 0, 0, 0, t5, -t17, t6, -g(1) * t51 - g(2) * (t48 + t66) - (t74 * t42 + t65) * t86, 0, 0, 0, 0, 0, 0, g(1) * t58 - g(2) * t4, -g(1) * t57 - g(2) * t3, t17, -g(1) * (-t13 * pkin(5) + t51) - g(2) * (pkin(8) * t79 + t73 * t15 + t48) - ((-pkin(8) + t74) * t42 + t65) * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, -t10, -g(1) * (-pkin(2) * t79 + t24) - g(2) * (-pkin(2) * t80 + t21) - g(3) * t70, 0, 0, 0, 0, 0, 0, -t10, -t8, t7, -g(1) * t24 - g(2) * t21 - g(3) * t64 + (pkin(3) * t40 + pkin(2) + t67) * t87, 0, 0, 0, 0, 0, 0, t7, t10, t8, -g(1) * t71 - g(2) * t72 - g(3) * t54 + (-t59 + t67) * t87, 0, 0, 0, 0, 0, 0, t9 * t56, t45 * t52 - t55 * t87, -t10, -g(1) * (pkin(8) * t77 + t71) - g(2) * (pkin(8) * t78 + t72) - g(3) * (t45 * t39 * pkin(5) + t54) + (-g(3) * pkin(8) + t18 * (t73 * t39 - t59)) * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 + g(2) * t57 + t42 * t52, g(1) * t4 + g(2) * t58 + t56 * t83, 0, 0;];
taug_reg  = t11;
