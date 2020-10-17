% Calculate minimal parameter regressor of gravitation load for
% S6RRRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPPR10_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR10_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 07:05:29
% EndTime: 2019-05-07 07:05:31
% DurationCPUTime: 0.56s
% Computational Cost: add. (428->125), mult. (1026->194), div. (0->0), fcn. (1267->12), ass. (0->62)
t53 = sin(qJ(3));
t56 = cos(qJ(3));
t96 = -pkin(3) * t56 - qJ(4) * t53 - pkin(2);
t55 = sin(qJ(1));
t54 = sin(qJ(2));
t77 = cos(pkin(6));
t69 = t54 * t77;
t89 = cos(qJ(2));
t90 = cos(qJ(1));
t31 = t55 * t89 + t90 * t69;
t51 = sin(pkin(6));
t74 = t51 * t90;
t13 = t31 * t53 + t56 * t74;
t62 = t77 * t89;
t30 = t55 * t54 - t90 * t62;
t49 = pkin(11) + qJ(6);
t46 = sin(t49);
t47 = cos(t49);
t95 = t13 * t46 + t30 * t47;
t94 = t13 * t47 - t30 * t46;
t93 = pkin(4) + pkin(9);
t91 = g(3) * t51;
t86 = t46 * t53;
t85 = t47 * t53;
t50 = sin(pkin(11));
t84 = t50 * t53;
t83 = t51 * t54;
t82 = t51 * t55;
t81 = t51 * t56;
t52 = cos(pkin(11));
t80 = t52 * t53;
t78 = qJ(5) * t56;
t76 = t96 * t30;
t32 = t90 * t54 + t55 * t62;
t75 = t96 * t32;
t73 = t51 * t89;
t72 = t53 * t89;
t71 = t56 * t89;
t14 = t31 * t56 - t53 * t74;
t70 = -t13 * pkin(3) + t14 * qJ(4);
t33 = -t55 * t69 + t90 * t89;
t17 = t33 * t53 - t55 * t81;
t18 = t33 * t56 + t53 * t82;
t68 = -t17 * pkin(3) + t18 * qJ(4);
t28 = t53 * t83 - t77 * t56;
t29 = t77 * t53 + t54 * t81;
t67 = -t28 * pkin(3) + t29 * qJ(4);
t66 = pkin(2) * t73 + pkin(9) * t83 + (pkin(3) * t71 + qJ(4) * t72) * t51;
t65 = -g(1) * t13 + g(2) * t17;
t64 = -g(1) * t14 + g(2) * t18;
t63 = g(1) * t30 - g(2) * t32;
t61 = t90 * pkin(1) + t33 * pkin(2) + t18 * pkin(3) + pkin(8) * t82 + t17 * qJ(4);
t2 = g(1) * t17 + g(2) * t13 + g(3) * t28;
t60 = g(1) * t18 + g(2) * t14 + g(3) * t29;
t59 = -t55 * pkin(1) - t31 * pkin(2) - pkin(3) * t14 + pkin(8) * t74 - qJ(4) * t13;
t58 = g(1) * t33 + g(2) * t31 + g(3) * t83;
t57 = -g(1) * t32 - g(2) * t30 + g(3) * t73;
t8 = t57 * t56;
t7 = t57 * t53;
t6 = t17 * t46 + t32 * t47;
t5 = t17 * t47 - t32 * t46;
t1 = [0, g(1) * t55 - g(2) * t90, g(1) * t90 + g(2) * t55, 0, 0, 0, 0, 0, g(1) * t31 - g(2) * t33, -t63, 0, 0, 0, 0, 0, -t64, t65, t63, t64, -t65, -g(1) * (-t30 * pkin(9) + t59) - g(2) * (t32 * pkin(9) + t61) -g(1) * (-t13 * t50 - t30 * t52) - g(2) * (t17 * t50 + t32 * t52) -g(1) * (-t13 * t52 + t30 * t50) - g(2) * (t17 * t52 - t32 * t50) -t64, -g(1) * (-qJ(5) * t14 - t93 * t30 + t59) - g(2) * (t18 * qJ(5) + t93 * t32 + t61) 0, 0, 0, 0, 0, g(1) * t95 - g(2) * t6, g(1) * t94 - g(2) * t5; 0, 0, 0, 0, 0, 0, 0, 0, -t57, t58, 0, 0, 0, 0, 0, -t8, t7, -t58, t8, -t7, -g(1) * (t33 * pkin(9) + t75) - g(2) * (t31 * pkin(9) + t76) - g(3) * t66, -g(1) * (-t32 * t84 + t33 * t52) - g(2) * (-t30 * t84 + t31 * t52) - (t50 * t72 + t52 * t54) * t91, -g(1) * (-t32 * t80 - t33 * t50) - g(2) * (-t30 * t80 - t31 * t50) - (-t50 * t54 + t52 * t72) * t91, -t8, -g(1) * (-t32 * t78 + t93 * t33 + t75) - g(2) * (-t30 * t78 + t93 * t31 + t76) - g(3) * ((pkin(4) * t54 + qJ(5) * t71) * t51 + t66) 0, 0, 0, 0, 0, -g(1) * (-t32 * t86 + t33 * t47) - g(2) * (-t30 * t86 + t31 * t47) - (t46 * t72 + t47 * t54) * t91, -g(1) * (-t32 * t85 - t33 * t46) - g(2) * (-t30 * t85 - t31 * t46) - (-t46 * t54 + t47 * t72) * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t60, 0, -t2, -t60, -g(1) * t68 - g(2) * t70 - g(3) * t67, -t60 * t50, -t60 * t52, t2, -g(1) * (-t17 * qJ(5) + t68) - g(2) * (-t13 * qJ(5) + t70) - g(3) * (-t28 * qJ(5) + t67) 0, 0, 0, 0, 0, -t60 * t46, -t60 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 - g(2) * t94 - g(3) * (t28 * t47 + t46 * t73) g(1) * t6 + g(2) * t95 - g(3) * (-t28 * t46 + t47 * t73);];
taug_reg  = t1;
