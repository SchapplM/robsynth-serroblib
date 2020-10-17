% Calculate minimal parameter regressor of gravitation load for
% S6RRRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [6x33]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRP7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP7_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP7_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 05:31:35
% EndTime: 2019-05-08 05:31:37
% DurationCPUTime: 0.49s
% Computational Cost: add. (504->100), mult. (865->164), div. (0->0), fcn. (1057->12), ass. (0->64)
t44 = sin(qJ(2));
t45 = sin(qJ(1));
t48 = cos(qJ(2));
t49 = cos(qJ(1));
t63 = cos(pkin(6));
t58 = t49 * t63;
t26 = t44 * t58 + t45 * t48;
t39 = qJ(3) + qJ(4);
t37 = sin(t39);
t38 = cos(t39);
t40 = sin(pkin(6));
t66 = t40 * t49;
t14 = t26 * t38 - t37 * t66;
t25 = t44 * t45 - t48 * t58;
t42 = sin(qJ(5));
t46 = cos(qJ(5));
t86 = t14 * t42 - t25 * t46;
t74 = t25 * t42;
t85 = t14 * t46 + t74;
t59 = t45 * t63;
t28 = -t44 * t59 + t48 * t49;
t43 = sin(qJ(3));
t47 = cos(qJ(3));
t67 = t40 * t47;
t19 = -t28 * t43 + t45 * t67;
t56 = t26 * t43 + t47 * t66;
t69 = t40 * t44;
t84 = g(2) * t56 - g(3) * (-t43 * t69 + t47 * t63) - g(1) * t19;
t24 = t37 * t63 + t38 * t69;
t64 = t46 * t48;
t68 = t40 * t45;
t18 = t28 * t38 + t37 * t68;
t27 = t44 * t49 + t48 * t59;
t7 = -t18 * t42 + t27 * t46;
t83 = g(2) * t86 - g(3) * (-t24 * t42 - t40 * t64) - g(1) * t7;
t82 = g(1) * t49 + g(2) * t45;
t78 = g(2) * t25;
t77 = g(2) * t26;
t75 = g(3) * t40;
t72 = t27 * t42;
t71 = t38 * t42;
t70 = t38 * t46;
t65 = t42 * t48;
t50 = -pkin(10) - pkin(9);
t61 = pkin(5) * t42 - t50;
t60 = t26 * t47 - t43 * t66;
t13 = t26 * t37 + t38 * t66;
t17 = t28 * t37 - t38 * t68;
t57 = -g(1) * t13 + g(2) * t17;
t35 = pkin(5) * t46 + pkin(4);
t36 = pkin(3) * t47 + pkin(2);
t41 = -qJ(6) - pkin(11);
t55 = t35 * t38 - t37 * t41 + t36;
t23 = t37 * t69 - t38 * t63;
t3 = g(1) * t17 + g(2) * t13 + g(3) * t23;
t5 = g(1) * t18 + g(2) * t14 + g(3) * t24;
t54 = -g(1) * t27 + t48 * t75 - t78;
t53 = -g(1) * (-t17 * t35 - t18 * t41) - g(2) * (-t13 * t35 - t14 * t41) - g(3) * (-t23 * t35 - t24 * t41);
t20 = t28 * t47 + t43 * t68;
t8 = t18 * t46 + t72;
t6 = t54 * t37;
t2 = t3 * t46;
t1 = t3 * t42;
t4 = [0, g(1) * t45 - g(2) * t49, t82, 0, 0, 0, 0, 0, g(1) * t26 - g(2) * t28, -g(1) * t25 + g(2) * t27, 0, 0, 0, 0, 0, g(1) * t60 - g(2) * t20, -g(1) * t56 - g(2) * t19, 0, 0, 0, 0, 0, g(1) * t14 - g(2) * t18, t57, 0, 0, 0, 0, 0, g(1) * t85 - g(2) * t8, -g(1) * t86 - g(2) * t7, -t57, -g(1) * (-t45 * pkin(1) - pkin(5) * t74 + t13 * t41 - t14 * t35 + t25 * t50 - t26 * t36) - g(2) * (t49 * pkin(1) + pkin(5) * t72 - t17 * t41 + t18 * t35 - t27 * t50 + t28 * t36) - t82 * t40 * (pkin(3) * t43 + pkin(8)); 0, 0, 0, 0, 0, 0, 0, 0, -t54, g(1) * t28 + g(3) * t69 + t77, 0, 0, 0, 0, 0, -t54 * t47, t54 * t43, 0, 0, 0, 0, 0, -t54 * t38, t6, 0, 0, 0, 0, 0, -g(1) * (-t27 * t70 + t28 * t42) - g(2) * (-t25 * t70 + t26 * t42) - (t38 * t64 + t42 * t44) * t75, -g(1) * (t27 * t71 + t28 * t46) - g(2) * (t25 * t71 + t26 * t46) - (-t38 * t65 + t44 * t46) * t75, -t6, -g(1) * (-t27 * t55 + t28 * t61) - t61 * t77 + t55 * t78 - (t44 * t61 + t48 * t55) * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, g(1) * t20 + g(2) * t60 - g(3) * (-t43 * t63 - t44 * t67) 0, 0, 0, 0, 0, t3, t5, 0, 0, 0, 0, 0, t2, -t1, -t5, pkin(3) * t84 + t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t5, 0, 0, 0, 0, 0, t2, -t1, -t5, t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, g(1) * t8 + g(2) * t85 - g(3) * (-t24 * t46 + t40 * t65) 0, t83 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3;];
taug_reg  = t4;
