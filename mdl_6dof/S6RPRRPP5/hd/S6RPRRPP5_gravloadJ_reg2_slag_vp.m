% Calculate inertial parameters regressor of gravitation load for
% S6RPRRPP5
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
% Datum: 2019-03-09 04:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRPP5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:38:25
% EndTime: 2019-05-05 21:38:26
% DurationCPUTime: 0.36s
% Computational Cost: add. (368->92), mult. (498->116), div. (0->0), fcn. (521->8), ass. (0->57)
t43 = sin(qJ(1));
t45 = cos(qJ(1));
t23 = g(1) * t45 + g(2) * t43;
t38 = pkin(9) + qJ(3);
t35 = sin(t38);
t77 = t23 * t35;
t36 = cos(t38);
t61 = t36 * pkin(3) + t35 * pkin(8);
t76 = -pkin(4) - pkin(5);
t41 = -pkin(7) - qJ(2);
t74 = g(2) * t41;
t42 = sin(qJ(4));
t72 = t35 * t42;
t71 = t35 * t43;
t44 = cos(qJ(4));
t70 = t35 * t44;
t69 = t35 * t45;
t68 = t36 * t44;
t67 = t36 * t45;
t66 = t43 * t42;
t65 = t43 * t44;
t64 = t45 * t41;
t63 = t45 * t42;
t62 = t45 * t44;
t60 = qJ(5) * t42;
t59 = qJ(6) * t36;
t58 = t35 * qJ(6);
t40 = cos(pkin(9));
t32 = t40 * pkin(2) + pkin(1);
t21 = t45 * t32;
t57 = pkin(3) * t67 + pkin(8) * t69 + t21;
t56 = -pkin(3) - t60;
t15 = t36 * t66 + t62;
t16 = t36 * t65 - t63;
t55 = -t15 * pkin(4) + t16 * qJ(5);
t17 = t36 * t63 - t65;
t18 = t36 * t62 + t66;
t54 = -t17 * pkin(4) + t18 * qJ(5);
t53 = pkin(4) * t68 + t36 * t60 + t61;
t4 = g(1) * t15 - g(2) * t17;
t22 = g(1) * t43 - g(2) * t45;
t52 = -t32 - t61;
t50 = -t16 * pkin(4) - t15 * qJ(5) - t64;
t49 = t18 * pkin(4) + t17 * qJ(5) + t57;
t2 = g(1) * t17 + g(2) * t15 + g(3) * t72;
t47 = g(1) * t18 + g(2) * t16 + g(3) * t70;
t8 = -g(3) * t36 + t77;
t46 = (-g(1) * t52 + t74) * t43;
t27 = pkin(8) * t67;
t24 = t43 * t36 * pkin(8);
t19 = qJ(5) * t70;
t14 = g(1) * t71 - g(2) * t69;
t9 = g(3) * t35 + t23 * t36;
t7 = t8 * t44;
t6 = t8 * t42;
t5 = g(1) * t16 - g(2) * t18;
t1 = [0, 0, 0, 0, 0, 0, t22, t23, 0, 0, 0, 0, 0, 0, 0, 0, t22 * t40, -t22 * sin(pkin(9)) -t23, -g(1) * (-t43 * pkin(1) + t45 * qJ(2)) - g(2) * (t45 * pkin(1) + t43 * qJ(2)) 0, 0, 0, 0, 0, 0, t22 * t36, -t14, -t23, -g(1) * (-t43 * t32 - t64) - g(2) * (-t43 * t41 + t21) 0, 0, 0, 0, 0, 0, t5, -t4, t14, g(1) * t64 - g(2) * t57 + t46, 0, 0, 0, 0, 0, 0, t5, t14, t4, -g(1) * t50 - g(2) * t49 + t46, 0, 0, 0, 0, 0, 0, t5, t4, -t14, -g(1) * (-t16 * pkin(5) + t50) - g(2) * (t18 * pkin(5) - t45 * t58 + t49) + (-g(1) * (t52 + t58) + t74) * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t9, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t6, -t9, -g(1) * (-pkin(3) * t69 + t27) - g(2) * (-pkin(3) * t71 + t24) - g(3) * t61, 0, 0, 0, 0, 0, 0, t7, -t9, t6, -g(1) * t27 - g(2) * t24 - g(3) * t53 + (pkin(4) * t44 - t56) * t77, 0, 0, 0, 0, 0, 0, t7, t6, t9, -g(1) * (-t45 * t59 + t27) - g(2) * (-t43 * t59 + t24) - g(3) * (pkin(5) * t68 + t53) + (g(3) * qJ(6) + t23 * (-t76 * t44 - t56)) * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t47, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, -t47, -g(1) * t54 - g(2) * t55 - g(3) * (-pkin(4) * t72 + t19) 0, 0, 0, 0, 0, 0, t2, -t47, 0, -g(1) * (-t17 * pkin(5) + t54) - g(2) * (-t15 * pkin(5) + t55) - g(3) * (t76 * t72 + t19); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8;];
taug_reg  = t1;
