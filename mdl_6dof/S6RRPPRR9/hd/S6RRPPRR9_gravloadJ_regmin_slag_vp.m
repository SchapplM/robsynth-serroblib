% Calculate minimal parameter regressor of gravitation load for
% S6RRPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPPRR9_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR9_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t35 = sin(qJ(2));
t36 = sin(qJ(1));
t39 = cos(qJ(2));
t40 = cos(qJ(1));
t51 = cos(pkin(6));
t49 = t40 * t51;
t20 = t36 * t35 - t39 * t49;
t33 = sin(qJ(6));
t37 = cos(qJ(6));
t21 = t35 * t49 + t36 * t39;
t34 = sin(qJ(5));
t38 = cos(qJ(5));
t32 = sin(pkin(6));
t57 = t32 * t40;
t44 = -t21 * t34 + t38 * t57;
t66 = -t20 * t37 + t33 * t44;
t65 = t20 * t33 + t37 * t44;
t64 = g(3) * t32;
t61 = t32 * t35;
t60 = t32 * t36;
t59 = t32 * t38;
t58 = t32 * t39;
t56 = t33 * t34;
t55 = t33 * t39;
t54 = t34 * t37;
t53 = t37 * t39;
t52 = pkin(2) * t58 + qJ(3) * t61;
t50 = t36 * t51;
t48 = -t20 * pkin(2) + t21 * qJ(3);
t22 = t40 * t35 + t39 * t50;
t23 = -t35 * t50 + t40 * t39;
t47 = -t22 * pkin(2) + t23 * qJ(3);
t6 = g(1) * t20 - g(2) * t22;
t7 = g(1) * t21 - g(2) * t23;
t46 = g(1) * t40 + g(2) * t36;
t45 = t40 * pkin(1) + t23 * pkin(2) + pkin(8) * t60 + t22 * qJ(3);
t10 = t21 * t38 + t34 * t57;
t8 = t23 * t38 - t34 * t60;
t43 = g(1) * t8 + g(2) * t10 + g(3) * (-t51 * t34 + t35 * t59);
t42 = -t36 * pkin(1) - t21 * pkin(2) + pkin(8) * t57 - t20 * qJ(3);
t3 = -g(1) * t22 - g(2) * t20 + g(3) * t58;
t41 = g(1) * t23 + g(2) * t21 + g(3) * t61;
t24 = t46 * t32;
t19 = t34 * t61 + t51 * t38;
t9 = t23 * t34 + t36 * t59;
t2 = -t22 * t33 + t9 * t37;
t1 = -t22 * t37 - t9 * t33;
t4 = [0, g(1) * t36 - g(2) * t40, t46, 0, 0, 0, 0, 0, t7, -t6, -t24, -t7, t6, -g(1) * t42 - g(2) * t45, -t24, t6, t7, -g(1) * (pkin(3) * t57 - t21 * qJ(4) + t42) - g(2) * (pkin(3) * t60 + t23 * qJ(4) + t45) 0, 0, 0, 0, 0, -g(1) * t44 - g(2) * t9, g(1) * t10 - g(2) * t8, 0, 0, 0, 0, 0, -g(1) * t65 - g(2) * t2, g(1) * t66 - g(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, -t3, t41, 0, t3, -t41, -g(1) * t47 - g(2) * t48 - g(3) * t52, 0, -t41, -t3, -g(1) * (-t22 * qJ(4) + t47) - g(2) * (-t20 * qJ(4) + t48) - g(3) * (qJ(4) * t58 + t52) 0, 0, 0, 0, 0, -t3 * t34, -t3 * t38, 0, 0, 0, 0, 0, -g(1) * (-t22 * t54 - t23 * t33) - g(2) * (-t20 * t54 - t21 * t33) - (-t33 * t35 + t34 * t53) * t64, -g(1) * (t22 * t56 - t23 * t37) - g(2) * (t20 * t56 - t21 * t37) - (-t34 * t55 - t35 * t37) * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, g(1) * t9 - g(2) * t44 + g(3) * t19, 0, 0, 0, 0, 0, -t43 * t37, t43 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 - g(2) * t66 - g(3) * (-t19 * t33 + t32 * t53) g(1) * t2 - g(2) * t65 - g(3) * (-t19 * t37 - t32 * t55);];
taug_reg  = t4;
