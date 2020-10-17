% Calculate minimal parameter regressor of gravitation load for
% S6PRRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% taug_reg [6x29]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRRP4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 09:55:42
% EndTime: 2019-05-05 09:55:43
% DurationCPUTime: 0.41s
% Computational Cost: add. (553->100), mult. (1101->159), div. (0->0), fcn. (1387->12), ass. (0->63)
t46 = sin(pkin(11));
t50 = sin(qJ(2));
t53 = cos(qJ(2));
t67 = cos(pkin(11));
t68 = cos(pkin(6));
t63 = t68 * t67;
t31 = t46 * t50 - t53 * t63;
t65 = t46 * t68;
t33 = t67 * t50 + t53 * t65;
t85 = -g(1) * t33 - g(2) * t31;
t32 = t46 * t53 + t50 * t63;
t34 = -t50 * t65 + t67 * t53;
t49 = sin(qJ(3));
t52 = cos(qJ(3));
t47 = sin(pkin(6));
t64 = t47 * t67;
t73 = t47 * t52;
t74 = t47 * t50;
t59 = g(3) * (-t49 * t74 + t68 * t52) + g(2) * (-t32 * t49 - t52 * t64) + g(1) * (-t34 * t49 + t46 * t73);
t79 = g(3) * t47;
t48 = sin(qJ(4));
t78 = t32 * t48;
t77 = t34 * t48;
t45 = qJ(4) + qJ(5);
t43 = sin(t45);
t76 = t43 * t52;
t44 = cos(t45);
t75 = t44 * t52;
t72 = t47 * t53;
t71 = t48 * t52;
t51 = cos(qJ(4));
t70 = t51 * t52;
t69 = t52 * t53;
t66 = t43 * t72;
t24 = t46 * t47 * t49 + t34 * t52;
t11 = t24 * t43 - t33 * t44;
t36 = t68 * t49 + t50 * t73;
t19 = t36 * t43 + t44 * t72;
t22 = t32 * t52 - t49 * t64;
t9 = t22 * t43 - t31 * t44;
t1 = g(1) * t11 + g(2) * t9 + g(3) * t19;
t10 = t22 * t44 + t31 * t43;
t12 = t24 * t44 + t33 * t43;
t20 = t36 * t44 - t66;
t3 = g(1) * t12 + g(2) * t10 + g(3) * t20;
t14 = -t31 * t76 - t32 * t44;
t16 = -t33 * t76 - t34 * t44;
t25 = -t44 * t74 + t52 * t66;
t60 = g(1) * t16 + g(2) * t14 + g(3) * t25;
t58 = g(1) * t24 + g(2) * t22 + g(3) * t36;
t57 = g(3) * t72 + t85;
t56 = -g(1) * (-t11 * pkin(5) + t12 * qJ(6)) - g(2) * (-t9 * pkin(5) + t10 * qJ(6)) - g(3) * (-t19 * pkin(5) + t20 * qJ(6));
t55 = -g(1) * (-t24 * t48 + t33 * t51) - g(2) * (-t22 * t48 + t31 * t51) - g(3) * (-t36 * t48 - t51 * t72);
t54 = -pkin(10) - pkin(9);
t42 = t51 * pkin(4) + pkin(3);
t26 = (t43 * t50 + t44 * t69) * t47;
t17 = -t33 * t75 + t34 * t43;
t15 = -t31 * t75 + t32 * t43;
t13 = t57 * t49;
t6 = t59 * t44;
t5 = t59 * t43;
t4 = -g(1) * t17 - g(2) * t15 - g(3) * t26;
t2 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, -t57, g(1) * t34 + g(2) * t32 + g(3) * t74, 0, 0, 0, 0, 0, -t57 * t52, t13, 0, 0, 0, 0, 0, -g(1) * (-t33 * t70 + t77) - g(2) * (-t31 * t70 + t78) - (t48 * t50 + t51 * t69) * t79, -g(1) * (t33 * t71 + t34 * t51) - g(2) * (t31 * t71 + t32 * t51) - (-t48 * t69 + t50 * t51) * t79, 0, 0, 0, 0, 0, t4, t60, t4, -t13, -t60, -g(1) * (pkin(4) * t77 + t17 * pkin(5) + t34 * pkin(8) + t16 * qJ(6)) - g(2) * (pkin(4) * t78 + t15 * pkin(5) + t32 * pkin(8) + t14 * qJ(6)) - g(3) * (t26 * pkin(5) + t25 * qJ(6)) - (pkin(4) * t48 + pkin(8)) * t50 * t79 + (-t53 * t79 - t85) * (t42 * t52 - t49 * t54 + pkin(2)); 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, t58, 0, 0, 0, 0, 0, -t59 * t51, t59 * t48, 0, 0, 0, 0, 0, -t6, t5, -t6, -t58, -t5, t58 * t54 - t59 * (pkin(5) * t44 + qJ(6) * t43 + t42); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, -g(1) * (-t24 * t51 - t33 * t48) - g(2) * (-t22 * t51 - t31 * t48) - g(3) * (-t36 * t51 + t48 * t72) 0, 0, 0, 0, 0, t1, t3, t1, 0, -t3, t55 * pkin(4) + t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t3, t1, 0, -t3, t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t2;
