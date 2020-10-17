% Calculate inertial parameters regressor of gravitation load for
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
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRPP3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:27:23
% EndTime: 2019-05-05 21:27:24
% DurationCPUTime: 0.34s
% Computational Cost: add. (394->92), mult. (485->111), div. (0->0), fcn. (510->8), ass. (0->54)
t36 = qJ(1) + pkin(9);
t31 = sin(t36);
t32 = cos(t36);
t18 = g(1) * t32 + g(2) * t31;
t38 = sin(qJ(3));
t73 = t18 * t38;
t72 = g(1) * t31;
t33 = t38 * pkin(8);
t41 = cos(qJ(3));
t34 = t41 * pkin(3);
t69 = t31 * t41;
t68 = t32 * t38;
t67 = t32 * t41;
t37 = sin(qJ(4));
t66 = t37 * t38;
t65 = t37 * t41;
t40 = cos(qJ(4));
t64 = t38 * t40;
t63 = t40 * t41;
t62 = -pkin(4) - qJ(6);
t61 = t34 + t33;
t60 = qJ(5) * t37;
t42 = cos(qJ(1));
t59 = t42 * pkin(1) + t32 * pkin(2) + t31 * pkin(7);
t58 = -pkin(2) - t34;
t39 = sin(qJ(1));
t57 = -t39 * pkin(1) + t32 * pkin(7);
t56 = -pkin(3) - t60;
t13 = t31 * t65 + t32 * t40;
t14 = t31 * t63 - t32 * t37;
t55 = -t13 * pkin(4) + t14 * qJ(5);
t15 = -t31 * t40 + t32 * t65;
t16 = t31 * t37 + t32 * t63;
t54 = -t15 * pkin(4) + t16 * qJ(5);
t53 = pkin(4) * t63 + t41 * t60 + t61;
t52 = pkin(3) * t67 + pkin(8) * t68 + t59;
t4 = g(1) * t13 - g(2) * t15;
t5 = g(1) * t14 - g(2) * t16;
t51 = -g(2) * t32 + t72;
t50 = g(1) * t39 - g(2) * t42;
t48 = -t14 * pkin(4) - t13 * qJ(5) + t57;
t46 = (t58 - t33) * t72;
t2 = g(1) * t15 + g(2) * t13 + g(3) * t66;
t45 = g(1) * t16 + g(2) * t14 + g(3) * t64;
t44 = t16 * pkin(4) + t15 * qJ(5) + t52;
t43 = -g(3) * t41 + t73;
t25 = qJ(5) * t64;
t21 = pkin(8) * t67;
t19 = pkin(8) * t69;
t17 = t51 * t38;
t8 = g(3) * t38 + t18 * t41;
t7 = t43 * t40;
t6 = t43 * t37;
t1 = [0, 0, 0, 0, 0, 0, t50, g(1) * t42 + g(2) * t39, 0, 0, 0, 0, 0, 0, 0, 0, t51, t18, 0, t50 * pkin(1), 0, 0, 0, 0, 0, 0, t51 * t41, -t17, -t18, -g(1) * (-t31 * pkin(2) + t57) - g(2) * t59, 0, 0, 0, 0, 0, 0, t5, -t4, t17, -g(1) * t57 - g(2) * t52 - t46, 0, 0, 0, 0, 0, 0, t17, -t5, t4, -g(1) * t48 - g(2) * t44 - t46, 0, 0, 0, 0, 0, 0, t17, t4, t5, -g(1) * (-t14 * qJ(6) + t48) - g(2) * (pkin(5) * t68 + t16 * qJ(6) + t44) - ((-pkin(5) - pkin(8)) * t38 + t58) * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t8, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t6, -t8, -g(1) * (-pkin(3) * t68 + t21) - g(2) * (-t31 * t38 * pkin(3) + t19) - g(3) * t61, 0, 0, 0, 0, 0, 0, -t8, -t7, t6, -g(1) * t21 - g(2) * t19 - g(3) * t53 + (pkin(4) * t40 - t56) * t73, 0, 0, 0, 0, 0, 0, -t8, t6, t7, -g(1) * (pkin(5) * t67 + t21) - g(2) * (pkin(5) * t69 + t19) - g(3) * (qJ(6) * t63 + t53) + (-g(3) * pkin(5) + t18 * (-t62 * t40 - t56)) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t45, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, -t45, -g(1) * t54 - g(2) * t55 - g(3) * (-pkin(4) * t66 + t25) 0, 0, 0, 0, 0, 0, 0, -t45, t2, -g(1) * (-t15 * qJ(6) + t54) - g(2) * (-t13 * qJ(6) + t55) - g(3) * (t62 * t66 + t25); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45;];
taug_reg  = t1;
