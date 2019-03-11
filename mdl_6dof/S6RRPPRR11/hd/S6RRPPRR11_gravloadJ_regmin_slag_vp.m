% Calculate minimal parameter regressor of gravitation load for
% S6RRPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPPRR11_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR11_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t37 = sin(qJ(2));
t38 = sin(qJ(1));
t40 = cos(qJ(2));
t41 = cos(qJ(1));
t53 = cos(pkin(6));
t51 = t41 * t53;
t20 = t37 * t51 + t38 * t40;
t36 = sin(qJ(6));
t39 = cos(qJ(6));
t19 = t38 * t37 - t40 * t51;
t32 = pkin(11) + qJ(5);
t29 = sin(t32);
t30 = cos(t32);
t34 = sin(pkin(6));
t57 = t34 * t41;
t45 = -t19 * t29 + t30 * t57;
t67 = t20 * t39 + t36 * t45;
t66 = -t20 * t36 + t39 * t45;
t65 = g(3) * t34;
t62 = t29 * t36;
t61 = t29 * t39;
t60 = t34 * t37;
t59 = t34 * t38;
t58 = t34 * t40;
t56 = t36 * t37;
t55 = t37 * t39;
t54 = pkin(2) * t58 + qJ(3) * t60;
t52 = t38 * t53;
t50 = -t19 * pkin(2) + t20 * qJ(3);
t21 = t41 * t37 + t40 * t52;
t22 = -t37 * t52 + t41 * t40;
t49 = -t21 * pkin(2) + t22 * qJ(3);
t48 = g(1) * t19 - g(2) * t21;
t6 = g(1) * t20 - g(2) * t22;
t47 = g(1) * t41 + g(2) * t38;
t46 = t41 * pkin(1) + t22 * pkin(2) + pkin(8) * t59 + t21 * qJ(3);
t9 = t19 * t30 + t29 * t57;
t7 = t21 * t30 - t29 * t59;
t44 = g(1) * t7 + g(2) * t9 + g(3) * (-t53 * t29 - t30 * t58);
t43 = -t38 * pkin(1) - t20 * pkin(2) + pkin(8) * t57 - t19 * qJ(3);
t3 = -g(1) * t21 - g(2) * t19 + g(3) * t58;
t42 = g(1) * t22 + g(2) * t20 + g(3) * t60;
t35 = cos(pkin(11));
t33 = sin(pkin(11));
t14 = -t29 * t58 + t53 * t30;
t8 = t21 * t29 + t30 * t59;
t2 = t22 * t36 + t8 * t39;
t1 = t22 * t39 - t8 * t36;
t4 = [0, g(1) * t38 - g(2) * t41, t47, 0, 0, 0, 0, 0, t6, -t48, -t47 * t34, -t6, t48, -g(1) * t43 - g(2) * t46, -g(1) * (-t19 * t33 + t35 * t57) - g(2) * (t21 * t33 + t35 * t59) -g(1) * (-t19 * t35 - t33 * t57) - g(2) * (t21 * t35 - t33 * t59) t6, -g(1) * (pkin(3) * t57 - t20 * qJ(4) + t43) - g(2) * (pkin(3) * t59 + t22 * qJ(4) + t46) 0, 0, 0, 0, 0, -g(1) * t45 - g(2) * t8, g(1) * t9 - g(2) * t7, 0, 0, 0, 0, 0, -g(1) * t66 - g(2) * t2, g(1) * t67 - g(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, -t3, t42, 0, t3, -t42, -g(1) * t49 - g(2) * t50 - g(3) * t54, -t42 * t33, -t42 * t35, -t3, -g(1) * (-t21 * qJ(4) + t49) - g(2) * (-t19 * qJ(4) + t50) - g(3) * (qJ(4) * t58 + t54) 0, 0, 0, 0, 0, -t42 * t29, -t42 * t30, 0, 0, 0, 0, 0, -g(1) * (-t21 * t36 + t22 * t61) - g(2) * (-t19 * t36 + t20 * t61) - (t29 * t55 + t36 * t40) * t65, -g(1) * (-t21 * t39 - t22 * t62) - g(2) * (-t19 * t39 - t20 * t62) - (-t29 * t56 + t39 * t40) * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, g(1) * t8 - g(2) * t45 + g(3) * t14, 0, 0, 0, 0, 0, -t44 * t39, t44 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 - g(2) * t67 - g(3) * (-t14 * t36 + t34 * t55) g(1) * t2 - g(2) * t66 - g(3) * (-t14 * t39 - t34 * t56);];
taug_reg  = t4;
