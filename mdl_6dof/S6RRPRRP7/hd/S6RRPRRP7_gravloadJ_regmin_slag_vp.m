% Calculate minimal parameter regressor of gravitation load for
% S6RRPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRP7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 18:16:18
% EndTime: 2019-05-06 18:16:19
% DurationCPUTime: 0.35s
% Computational Cost: add. (266->81), mult. (705->104), div. (0->0), fcn. (824->8), ass. (0->50)
t40 = sin(qJ(5));
t44 = cos(qJ(5));
t41 = sin(qJ(4));
t45 = cos(qJ(2));
t42 = sin(qJ(2));
t60 = cos(qJ(4));
t55 = t42 * t60;
t24 = -t45 * t41 + t55;
t43 = sin(qJ(1));
t17 = t24 * t43;
t46 = cos(qJ(1));
t58 = t45 * t46;
t19 = t41 * t58 - t46 * t55;
t23 = t42 * t41 + t45 * t60;
t48 = g(1) * t19 - g(2) * t17 + g(3) * t23;
t74 = t48 * (pkin(5) * t44 + qJ(6) * t40 + pkin(4));
t25 = g(1) * t46 + g(2) * t43;
t73 = t25 * t42;
t72 = t48 * t40;
t71 = t48 * t44;
t34 = t42 * qJ(3);
t57 = t45 * pkin(2) + t34;
t69 = pkin(2) * t42;
t67 = g(1) * t43;
t62 = g(3) * t24;
t61 = t45 * pkin(3);
t56 = qJ(3) * t45;
t54 = pkin(2) * t58 + t43 * pkin(7) + (pkin(1) + t34) * t46;
t20 = t23 * t46;
t13 = t20 * t40 + t43 * t44;
t18 = t23 * t43;
t9 = t18 * t40 - t46 * t44;
t53 = -g(1) * t9 + g(2) * t13;
t52 = -g(1) * t17 - g(2) * t19;
t51 = -g(2) * t46 + t67;
t10 = t18 * t44 + t46 * t40;
t50 = -pkin(1) - t57;
t8 = g(1) * t20 + g(2) * t18 + t62;
t1 = g(1) * t13 + g(2) * t9 + t40 * t62;
t14 = t20 * t44 - t43 * t40;
t47 = g(1) * t14 + g(2) * t10 + t44 * t62;
t37 = t46 * pkin(7);
t30 = t46 * t56;
t28 = t43 * t56;
t22 = t51 * t45;
t21 = t51 * t42;
t16 = g(3) * t42 + t25 * t45;
t15 = -g(3) * t45 + t73;
t6 = g(1) * t10 - g(2) * t14;
t2 = [0, t51, t25, 0, 0, 0, 0, 0, t22, -t21, t22, -t25, t21, -g(1) * t37 - g(2) * t54 - t50 * t67, 0, 0, 0, 0, 0, g(1) * t18 - g(2) * t20, -t52, 0, 0, 0, 0, 0, t6, t53, t6, t52, -t53, -g(1) * (-t18 * pkin(4) - pkin(5) * t10 - t46 * pkin(8) + t17 * pkin(9) - qJ(6) * t9 + t37) - g(2) * (pkin(3) * t58 + t20 * pkin(4) + t14 * pkin(5) + t19 * pkin(9) + t13 * qJ(6) + t54) + (-g(1) * (t50 - t61) + g(2) * pkin(8)) * t43; 0, 0, 0, 0, 0, 0, 0, 0, t15, t16, t15, 0, -t16, -g(1) * (-t46 * t69 + t30) - g(2) * (-t43 * t69 + t28) - g(3) * t57, 0, 0, 0, 0, 0, -t48, -t8, 0, 0, 0, 0, 0, -t71, t72, -t71, t8, -t72, -g(1) * (-t20 * pkin(9) + t30) - g(2) * (-t18 * pkin(9) + t28) - g(3) * (-t24 * pkin(9) + t57 + t61) + (pkin(2) + pkin(3)) * t73 - t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t8, 0, 0, 0, 0, 0, t71, -t72, t71, -t8, t72, -t8 * pkin(9) + t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t47, t1, 0, -t47, -g(1) * (-t13 * pkin(5) + t14 * qJ(6)) - g(2) * (-t9 * pkin(5) + t10 * qJ(6)) - (-pkin(5) * t40 + qJ(6) * t44) * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t2;
