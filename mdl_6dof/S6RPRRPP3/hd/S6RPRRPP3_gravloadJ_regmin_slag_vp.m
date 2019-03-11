% Calculate minimal parameter regressor of gravitation load for
% S6RPRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% 
% Output:
% taug_reg [6x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRPP3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t37 = sin(qJ(3));
t35 = qJ(1) + pkin(9);
t30 = sin(t35);
t31 = cos(t35);
t50 = g(1) * t31 + g(2) * t30;
t68 = t50 * t37;
t67 = g(1) * t30;
t32 = t37 * pkin(8);
t40 = cos(qJ(3));
t33 = t40 * pkin(3);
t64 = t30 * t40;
t63 = t31 * t37;
t62 = t31 * t40;
t36 = sin(qJ(4));
t61 = t36 * t37;
t60 = t36 * t40;
t39 = cos(qJ(4));
t59 = t37 * t39;
t58 = t39 * t40;
t57 = -pkin(4) - qJ(6);
t56 = qJ(5) * t36;
t55 = -pkin(2) - t33;
t54 = -pkin(3) - t56;
t13 = t30 * t60 + t31 * t39;
t14 = t30 * t58 - t31 * t36;
t53 = -t13 * pkin(4) + t14 * qJ(5);
t15 = -t30 * t39 + t31 * t60;
t16 = t30 * t36 + t31 * t58;
t52 = -t15 * pkin(4) + t16 * qJ(5);
t51 = pkin(4) * t58 + t40 * t56 + t32 + t33;
t4 = g(1) * t13 - g(2) * t15;
t5 = g(1) * t14 - g(2) * t16;
t49 = -g(2) * t31 + t67;
t38 = sin(qJ(1));
t41 = cos(qJ(1));
t48 = g(1) * t38 - g(2) * t41;
t46 = -t38 * pkin(1) - t14 * pkin(4) + t31 * pkin(7) - t13 * qJ(5);
t2 = g(1) * t15 + g(2) * t13 + g(3) * t61;
t44 = g(1) * t16 + g(2) * t14 + g(3) * t59;
t43 = t41 * pkin(1) + t31 * pkin(2) + pkin(3) * t62 + t16 * pkin(4) + t30 * pkin(7) + pkin(8) * t63 + t15 * qJ(5);
t42 = -g(3) * t40 + t68;
t24 = qJ(5) * t59;
t20 = pkin(8) * t62;
t18 = pkin(8) * t64;
t17 = t49 * t37;
t8 = g(3) * t37 + t50 * t40;
t7 = t42 * t39;
t6 = t42 * t36;
t1 = [0, t48, g(1) * t41 + g(2) * t38, t48 * pkin(1), 0, 0, 0, 0, 0, t49 * t40, -t17, 0, 0, 0, 0, 0, t5, -t4, t17, -t5, t4, -g(1) * t46 - g(2) * t43 - (t55 - t32) * t67, t17, t4, t5, -g(1) * (-t14 * qJ(6) + t46) - g(2) * (pkin(5) * t63 + t16 * qJ(6) + t43) - ((-pkin(5) - pkin(8)) * t37 + t55) * t67; 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, t8, 0, 0, 0, 0, 0, t7, -t6, -t8, -t7, t6, -g(1) * t20 - g(2) * t18 - g(3) * t51 + (pkin(4) * t39 - t54) * t68, -t8, t6, t7, -g(1) * (pkin(5) * t62 + t20) - g(2) * (pkin(5) * t64 + t18) - g(3) * (qJ(6) * t58 + t51) + (-g(3) * pkin(5) + t50 * (-t57 * t39 - t54)) * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t44, 0, -t2, -t44, -g(1) * t52 - g(2) * t53 - g(3) * (-pkin(4) * t61 + t24) 0, -t44, t2, -g(1) * (-t15 * qJ(6) + t52) - g(2) * (-t13 * qJ(6) + t53) - g(3) * (t57 * t61 + t24); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44;];
taug_reg  = t1;
