% Calculate minimal parameter regressor of gravitation load for
% S6RRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
% 
% Output:
% taug_reg [6x30]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPPR4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 04:46:09
% EndTime: 2019-05-07 04:46:10
% DurationCPUTime: 0.38s
% Computational Cost: add. (319->79), mult. (484->122), div. (0->0), fcn. (540->10), ass. (0->52)
t32 = qJ(3) + pkin(10);
t26 = sin(t32);
t27 = cos(t32);
t34 = sin(qJ(6));
t38 = cos(qJ(6));
t47 = t26 * t38 - t27 * t34;
t37 = sin(qJ(1));
t40 = cos(qJ(2));
t41 = cos(qJ(1));
t59 = t41 * t26;
t8 = -t37 * t27 + t40 * t59;
t58 = t41 * t27;
t9 = t37 * t26 + t40 * t58;
t50 = t9 * t34 - t8 * t38;
t36 = sin(qJ(2));
t66 = g(3) * t36;
t60 = t37 * t40;
t6 = t26 * t60 + t58;
t7 = t27 * t60 - t59;
t75 = -t7 * t34 + t6 * t38;
t80 = -g(1) * t50 + g(2) * t75 + t47 * t66;
t39 = cos(qJ(3));
t25 = t39 * pkin(3) + pkin(2);
t18 = t40 * t25;
t33 = -qJ(4) - pkin(8);
t79 = -t36 * t33 + t18;
t53 = g(1) * t41 + g(2) * t37;
t78 = -g(3) * t40 + t53 * t36;
t77 = -pkin(1) - t79;
t2 = t8 * t34 + t9 * t38;
t46 = t26 * t34 + t27 * t38;
t51 = t6 * t34 + t7 * t38;
t73 = g(1) * t2 + g(2) * t51 + t46 * t66;
t35 = sin(qJ(3));
t62 = t37 * t35;
t61 = t37 * t39;
t57 = t41 * t35;
t56 = t41 * t39;
t55 = t40 * t57;
t52 = g(1) * t37 - g(2) * t41;
t49 = pkin(4) * t27 + qJ(5) * t26;
t12 = t35 * t60 + t56;
t44 = pkin(3) * t62 + t37 * pkin(7) - t77 * t41;
t43 = g(1) * t8 + g(2) * t6 + t26 * t66;
t42 = pkin(3) * t57 + t41 * pkin(7) + t77 * t37;
t11 = t53 * t40 + t66;
t23 = pkin(3) * t61;
t16 = t52 * t36;
t15 = t40 * t56 + t62;
t14 = -t55 + t61;
t13 = -t39 * t60 + t57;
t1 = [0, t52, t53, 0, 0, 0, 0, 0, t52 * t40, -t16, 0, 0, 0, 0, 0, -g(1) * t13 - g(2) * t15, -g(1) * t12 - g(2) * t14, t16, -g(1) * t42 - g(2) * t44, g(1) * t7 - g(2) * t9, t16, g(1) * t6 - g(2) * t8, -g(1) * (-t7 * pkin(4) - t6 * qJ(5) + t42) - g(2) * (t9 * pkin(4) + t8 * qJ(5) + t44) 0, 0, 0, 0, 0, g(1) * t51 - g(2) * t2, g(1) * t75 + g(2) * t50; 0, 0, 0, 0, 0, 0, 0, 0, t78, t11, 0, 0, 0, 0, 0, t78 * t39, -t78 * t35, -t11, -g(3) * t79 + t53 * (t25 * t36 + t33 * t40) t78 * t27, -t11, t78 * t26, -g(3) * t18 + (-g(3) * t49 + t53 * t33) * t40 + (g(3) * t33 + t53 * (t25 + t49)) * t36, 0, 0, 0, 0, 0, t78 * t46, t78 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t14 + g(2) * t12 + t35 * t66, g(1) * t15 - g(2) * t13 + t39 * t66, 0, -g(1) * t23 + (g(2) * t56 + t11 * t35) * pkin(3), t43, 0, -g(1) * t9 - g(2) * t7 - t27 * t66, -g(1) * (-pkin(3) * t55 - t8 * pkin(4) + t9 * qJ(5) + t23) - g(2) * (-t12 * pkin(3) - t6 * pkin(4) + t7 * qJ(5)) - (-pkin(3) * t35 - pkin(4) * t26 + qJ(5) * t27) * t66, 0, 0, 0, 0, 0, t80, -t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, 0, 0, 0, -t78, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80, t73;];
taug_reg  = t1;
