% Calculate minimal parameter regressor of gravitation load for
% S6RRRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% 
% Output:
% taug_reg [6x30]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRP8_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP8_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 08:35:51
% EndTime: 2019-05-07 08:35:53
% DurationCPUTime: 0.53s
% Computational Cost: add. (241->94), mult. (605->133), div. (0->0), fcn. (680->8), ass. (0->65)
t42 = sin(qJ(1));
t44 = cos(qJ(3));
t45 = cos(qJ(2));
t40 = sin(qJ(3));
t46 = cos(qJ(1));
t66 = t46 * t40;
t19 = -t42 * t44 + t45 * t66;
t65 = t46 * t44;
t20 = t42 * t40 + t45 * t65;
t39 = sin(qJ(5));
t43 = cos(qJ(5));
t3 = t19 * t43 - t20 * t39;
t69 = t42 * t45;
t17 = t40 * t69 + t65;
t68 = t44 * t45;
t18 = t42 * t68 - t66;
t90 = t17 * t43 - t18 * t39;
t93 = g(1) * t3 + g(2) * t90;
t41 = sin(qJ(2));
t73 = t39 * t44;
t53 = -t40 * t43 + t73;
t91 = t53 * t41;
t92 = -g(3) * t91 + t93;
t56 = g(1) * t46 + g(2) * t42;
t88 = t41 * t56;
t8 = -g(3) * t45 + t88;
t79 = g(3) * t41;
t4 = t19 * t39 + t20 * t43;
t74 = t39 * t40;
t52 = t43 * t44 + t74;
t76 = t17 * t39;
t54 = t18 * t43 + t76;
t87 = g(1) * t4 + g(2) * t54 + t52 * t79;
t84 = g(1) * t42;
t33 = t41 * pkin(8);
t35 = t45 * pkin(2);
t32 = t43 * pkin(5) + pkin(4);
t77 = -pkin(3) - t32;
t72 = t40 * t41;
t71 = t41 * t44;
t70 = t41 * t46;
t67 = t45 * t46;
t64 = qJ(4) * t40;
t63 = -pkin(1) - t35;
t62 = pkin(5) * t39 + qJ(4);
t60 = pkin(3) * t68 + t45 * t64 + t33 + t35;
t59 = -t18 * pkin(3) + t46 * pkin(7) - t17 * qJ(4);
t58 = t46 * pkin(1) + pkin(2) * t67 + t20 * pkin(3) + t42 * pkin(7) + pkin(8) * t70;
t57 = g(1) * t17 - g(2) * t19;
t55 = -g(2) * t46 + t84;
t49 = g(3) * t53;
t2 = g(1) * t19 + g(2) * t17 + g(3) * t72;
t48 = g(1) * t20 + g(2) * t18 + g(3) * t71;
t38 = -qJ(6) - pkin(9);
t27 = pkin(8) * t67;
t24 = pkin(8) * t69;
t22 = qJ(4) * t71;
t21 = -g(2) * t70 + t41 * t84;
t15 = t19 * pkin(3);
t13 = t17 * pkin(3);
t9 = t45 * t56 + t79;
t7 = t8 * t44;
t6 = t8 * t40;
t5 = g(1) * t18 - g(2) * t20;
t1 = [0, t55, t56, 0, 0, 0, 0, 0, t55 * t45, -t21, 0, 0, 0, 0, 0, t5, -t57, t5, t21, t57, -g(1) * t59 - g(2) * (t19 * qJ(4) + t58) - (t63 - t33) * t84, 0, 0, 0, 0, 0, g(1) * t54 - g(2) * t4, g(1) * t90 - g(2) * t3, -t21, -g(1) * (-pkin(5) * t76 - t18 * t32 + t59) - g(2) * (t19 * t62 + t20 * t32 + t38 * t70 + t58) - ((-pkin(8) - t38) * t41 + t63) * t84; 0, 0, 0, 0, 0, 0, 0, 0, t8, t9, 0, 0, 0, 0, 0, t7, -t6, t7, -t9, t6, -g(1) * t27 - g(2) * t24 - g(3) * t60 + (pkin(3) * t44 + pkin(2) + t64) * t88, 0, 0, 0, 0, 0, t8 * t52, t45 * t49 - t56 * t91, t9, -g(1) * (t38 * t67 + t27) - g(2) * (t38 * t69 + t24) - g(3) * (pkin(5) * t45 * t74 + t32 * t68 + t60) + (-g(3) * t38 + t56 * (t40 * t62 - t44 * t77 + pkin(2))) * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t48, t2, 0, -t48, -g(1) * (t20 * qJ(4) - t15) - g(2) * (t18 * qJ(4) - t13) - g(3) * (-pkin(3) * t72 + t22) 0, 0, 0, 0, 0, t92, -t87, 0, -g(1) * (-t19 * t32 + t20 * t62 - t15) - g(2) * (-t17 * t32 + t18 * t62 - t13) - g(3) * t22 - (pkin(5) * t73 + t40 * t77) * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t92, t87, 0 (t41 * t49 - t93) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8;];
taug_reg  = t1;
