% Calculate minimal parameter regressor of gravitation load for
% S6RRPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRPR14_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR14_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 17:00:04
% EndTime: 2019-05-06 17:00:05
% DurationCPUTime: 0.35s
% Computational Cost: add. (279->96), mult. (729->146), div. (0->0), fcn. (891->10), ass. (0->52)
t39 = sin(qJ(2));
t40 = sin(qJ(1));
t43 = cos(qJ(2));
t44 = cos(qJ(1));
t61 = cos(pkin(6));
t58 = t44 * t61;
t25 = t39 * t58 + t40 * t43;
t59 = t40 * t61;
t27 = -t39 * t59 + t44 * t43;
t75 = -g(1) * t27 - g(2) * t25;
t37 = sin(qJ(6));
t41 = cos(qJ(6));
t24 = t40 * t39 - t43 * t58;
t38 = sin(qJ(4));
t42 = cos(qJ(4));
t36 = sin(pkin(6));
t64 = t36 * t44;
t50 = t24 * t42 + t38 * t64;
t74 = -t25 * t41 + t37 * t50;
t73 = t25 * t37 + t41 * t50;
t70 = g(3) * t36;
t67 = t36 * t39;
t66 = t36 * t40;
t65 = t36 * t43;
t63 = t37 * t42;
t62 = t41 * t42;
t60 = g(3) * (pkin(2) * t65 + qJ(3) * t67);
t26 = t44 * t39 + t43 * t59;
t9 = -t26 * t42 + t38 * t66;
t57 = g(1) * t50 + g(2) * t9;
t10 = t26 * t38 + t42 * t66;
t51 = -t24 * t38 + t42 * t64;
t56 = g(1) * t51 + g(2) * t10;
t55 = g(1) * t24 - g(2) * t26;
t8 = g(1) * t25 - g(2) * t27;
t54 = g(1) * t44 + g(2) * t40;
t53 = pkin(4) * t38 - qJ(5) * t42;
t52 = t44 * pkin(1) + t27 * pkin(2) + pkin(8) * t66 + t26 * qJ(3);
t22 = t61 * t38 + t42 * t65;
t48 = g(1) * t9 - g(2) * t50 + g(3) * t22;
t47 = -t40 * pkin(1) - t25 * pkin(2) + pkin(8) * t64 - t24 * qJ(3);
t23 = -t38 * t65 + t61 * t42;
t46 = g(1) * t10 - g(2) * t51 + g(3) * t23;
t6 = -g(1) * t26 - g(2) * t24 + g(3) * t65;
t45 = g(3) * t67 - t75;
t20 = t26 * pkin(2);
t18 = t24 * pkin(2);
t5 = t45 * t42;
t4 = t45 * t38;
t3 = t27 * t41 + t9 * t37;
t2 = -t27 * t37 + t9 * t41;
t1 = [0, g(1) * t40 - g(2) * t44, t54, 0, 0, 0, 0, 0, t8, -t55, -t54 * t36, -t8, t55, -g(1) * t47 - g(2) * t52, 0, 0, 0, 0, 0, -t56, t57, t8, t56, -t57, -g(1) * (pkin(3) * t64 + pkin(4) * t51 - t25 * pkin(9) + qJ(5) * t50 + t47) - g(2) * (pkin(3) * t66 + t10 * pkin(4) + t27 * pkin(9) + t9 * qJ(5) + t52) 0, 0, 0, 0, 0, -g(1) * t74 - g(2) * t3, -g(1) * t73 - g(2) * t2; 0, 0, 0, 0, 0, 0, 0, 0, -t6, t45, 0, t6, -t45, -g(1) * (t27 * qJ(3) - t20) - g(2) * (t25 * qJ(3) - t18) - t60, 0, 0, 0, 0, 0, -t4, -t5, -t6, t4, t5, -g(1) * (-t26 * pkin(9) - t20) - g(2) * (-t24 * pkin(9) - t18) - t60 - (pkin(9) * t43 + t53 * t39) * t70 + t75 * (qJ(3) + t53) 0, 0, 0, 0, 0, -g(1) * (-t26 * t41 - t27 * t63) - g(2) * (-t24 * t41 - t25 * t63) - (-t39 * t63 + t41 * t43) * t70, -g(1) * (t26 * t37 - t27 * t62) - g(2) * (t24 * t37 - t25 * t62) - (-t37 * t43 - t39 * t62) * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t46, 0, -t48, -t46, -g(1) * (-t9 * pkin(4) + t10 * qJ(5)) - g(2) * (pkin(4) * t50 - qJ(5) * t51) - g(3) * (-t22 * pkin(4) + t23 * qJ(5)) 0, 0, 0, 0, 0, -t46 * t37, -t46 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t2 + g(2) * t73 - g(3) * (t22 * t41 - t37 * t67) g(1) * t3 - g(2) * t74 - g(3) * (-t22 * t37 - t41 * t67);];
taug_reg  = t1;
