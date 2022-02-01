% Calculate inertial parameters regressor of gravitation load for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRR4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:16:59
% EndTime: 2022-01-23 09:16:59
% DurationCPUTime: 0.25s
% Computational Cost: add. (189->63), mult. (226->91), div. (0->0), fcn. (239->10), ass. (0->52)
t31 = pkin(9) + qJ(4);
t23 = cos(t31);
t33 = sin(pkin(8));
t35 = cos(pkin(8));
t40 = qJ(3) + pkin(6);
t34 = cos(pkin(9));
t55 = t34 * pkin(3) + pkin(2);
t61 = -(-pkin(7) - t40) * t33 + (pkin(4) * t23 + t55) * t35;
t22 = sin(t31);
t57 = g(3) * t33;
t37 = cos(qJ(1));
t43 = t37 * t23;
t36 = sin(qJ(1));
t50 = t36 * t22;
t7 = t35 * t50 + t43;
t44 = t37 * t22;
t49 = t36 * t23;
t9 = -t35 * t44 + t49;
t60 = -g(1) * t9 + g(2) * t7 + t22 * t57;
t32 = sin(pkin(9));
t56 = t32 * pkin(3);
t24 = qJ(5) + t31;
t20 = sin(t24);
t52 = t36 * t20;
t21 = cos(t24);
t51 = t36 * t21;
t48 = t36 * t32;
t47 = t36 * t34;
t46 = t37 * t20;
t45 = t37 * t21;
t42 = t37 * t32;
t41 = t37 * t34;
t25 = t36 * qJ(2);
t39 = t37 * pkin(1) + t25;
t18 = g(1) * t37 + g(2) * t36;
t17 = g(1) * t36 - g(2) * t37;
t27 = t37 * qJ(2);
t19 = qJ(2) + t56;
t16 = pkin(2) * t35 + t33 * qJ(3) + pkin(1);
t15 = pkin(4) * t22 + t56;
t13 = t17 * t33;
t12 = t40 * t33 + t55 * t35 + pkin(1);
t11 = g(3) * t35 - t18 * t33;
t10 = t35 * t43 + t50;
t8 = -t35 * t49 + t44;
t6 = t35 * t45 + t52;
t5 = -t35 * t46 + t51;
t4 = -t35 * t51 + t46;
t3 = t35 * t52 + t45;
t2 = g(1) * t6 - g(2) * t4 + t21 * t57;
t1 = -g(1) * t5 + g(2) * t3 + t20 * t57;
t14 = [0, 0, 0, 0, 0, 0, t17, t18, 0, 0, 0, 0, 0, 0, 0, 0, t17 * t35, -t13, -t18, -g(1) * (-t36 * pkin(1) + t27) - g(2) * t39, 0, 0, 0, 0, 0, 0, -g(1) * (-t35 * t47 + t42) - g(2) * (t35 * t41 + t48), -g(1) * (t35 * t48 + t41) - g(2) * (-t35 * t42 + t47), t13, -g(1) * (-t16 * t36 + t27) - g(2) * (t16 * t37 + t25), 0, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t10, -g(1) * t7 - g(2) * t9, t13, -g(1) * (-t12 * t36 + t19 * t37) - g(2) * (t12 * t37 + t19 * t36), 0, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5, t13, -g(1) * (t37 * t15 + t27) - g(2) * (t61 * t37 + t39) + (-g(1) * (-pkin(1) - t61) - g(2) * t15) * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, g(1) * t10 - g(2) * t8 + t23 * t57, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t60 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t14;
