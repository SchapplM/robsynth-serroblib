% Calculate minimal parameter regressor of gravitation load for
% S5RRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
% 
% Output:
% taug_reg [5x22]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPPP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:24:30
% EndTime: 2019-12-31 19:24:31
% DurationCPUTime: 0.31s
% Computational Cost: add. (225->80), mult. (636->116), div. (0->0), fcn. (727->8), ass. (0->65)
t46 = cos(pkin(5));
t85 = qJ(3) * t46 + pkin(7);
t47 = sin(qJ(2));
t44 = sin(pkin(5));
t68 = qJ(3) * t44;
t34 = t47 * t68;
t49 = cos(qJ(2));
t69 = t49 * pkin(2) + t34;
t43 = sin(pkin(8));
t78 = t47 * t43;
t64 = t46 * t78;
t45 = cos(pkin(8));
t72 = t49 * t45;
t23 = -t64 + t72;
t48 = sin(qJ(1));
t50 = cos(qJ(1));
t59 = g(1) * t50 + g(2) * t48;
t84 = g(3) * t47 + t49 * t59;
t83 = g(1) * t48;
t81 = t44 * t48;
t80 = t44 * t49;
t79 = t44 * t50;
t77 = t47 * t45;
t76 = t47 * t48;
t75 = t47 * t50;
t74 = t48 * t46;
t73 = t49 * t43;
t71 = t49 * t50;
t70 = t85 * t50;
t66 = pkin(2) * t76;
t65 = pkin(2) * t75;
t63 = t49 * t68;
t54 = -t46 * t72 + t78;
t17 = t54 * t48;
t53 = t46 * t73 + t77;
t18 = t53 * t48;
t28 = t48 * t63;
t62 = -t18 * pkin(3) - t17 * qJ(4) + t28;
t19 = t54 * t50;
t20 = t53 * t50;
t30 = t50 * t63;
t61 = -t20 * pkin(3) - t19 * qJ(4) + t30;
t60 = pkin(2) * t71 + (pkin(1) + t34) * t50 + t85 * t48;
t10 = t48 * t73 + (t47 * t74 + t79) * t45;
t22 = t46 * t77 + t73;
t12 = t22 * t50 - t45 * t81;
t2 = g(1) * t10 - g(2) * t12;
t11 = t23 * t48 - t43 * t79;
t13 = t43 * t81 + t45 * t71 - t50 * t64;
t3 = g(1) * t11 - g(2) * t13;
t58 = -g(2) * t50 + t83;
t57 = -t11 * pkin(3) - t10 * qJ(4) + t70;
t56 = t23 * pkin(3) + qJ(4) * t22 + t69;
t55 = -pkin(2) * t47 + pkin(4) * t80;
t4 = g(1) * t19 + g(2) * t17 - g(3) * t22;
t5 = g(1) * t20 + g(2) * t18 - g(3) * t23;
t52 = t13 * pkin(3) + qJ(4) * t12 + t60;
t51 = (-pkin(1) - t69) * t83;
t25 = t44 * t75 + t74;
t24 = t44 * t76 - t46 * t50;
t14 = t84 * t44;
t7 = g(1) * t24 - g(2) * t25;
t6 = -g(1) * t25 - g(2) * t24 + g(3) * t80;
t1 = -g(1) * t12 - g(2) * t10 - g(3) * t54;
t8 = [0, t58, t59, 0, 0, 0, 0, 0, t58 * t49, -t58 * t47, t3, -t2, t7, -g(1) * t70 - g(2) * t60 - t51, t7, -t3, t2, -g(1) * t57 - g(2) * t52 - t51, t7, t2, t3, -g(1) * (-t24 * pkin(4) - t11 * qJ(5) + t57) - g(2) * (pkin(4) * t25 + qJ(5) * t13 + t52) - t51; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t49 + t47 * t59, t84, t5, -t4, -t14, -g(1) * (t30 - t65) - g(2) * (t28 - t66) - g(3) * t69, -t14, -t5, t4, -g(1) * (t61 - t65) - g(2) * (t62 - t66) - g(3) * t56, -t14, t4, t5, -g(1) * (-t20 * qJ(5) + t50 * t55 + t61) - g(2) * (-t18 * qJ(5) + t48 * t55 + t62) - g(3) * (pkin(4) * t44 * t47 + qJ(5) * t23 + t56); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, t6, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t13 - g(2) * t11 - g(3) * t53;];
taug_reg = t8;
