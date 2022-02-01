% Calculate inertial parameters regressor of gravitation load for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:18
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:17:25
% EndTime: 2022-01-20 11:17:26
% DurationCPUTime: 0.27s
% Computational Cost: add. (298->63), mult. (283->86), div. (0->0), fcn. (294->10), ass. (0->53)
t37 = qJ(1) + qJ(2);
t32 = sin(t37);
t34 = cos(t37);
t42 = cos(qJ(4));
t39 = cos(pkin(9));
t40 = sin(qJ(4));
t54 = t39 * t40;
t11 = t32 * t54 + t34 * t42;
t13 = t32 * t42 - t34 * t54;
t38 = sin(pkin(9));
t62 = g(3) * t38;
t66 = -g(1) * t13 + g(2) * t11 + t40 * t62;
t64 = g(1) * t32;
t41 = sin(qJ(1));
t61 = t41 * pkin(1);
t60 = t32 * t39;
t59 = t32 * t40;
t58 = t34 * t38;
t57 = t34 * t39;
t56 = t34 * t40;
t55 = t38 * (-pkin(8) - pkin(7));
t53 = t39 * t42;
t52 = t34 * pkin(2) + t32 * qJ(3);
t27 = t34 * qJ(3);
t50 = -t32 * pkin(2) + t27;
t49 = pkin(3) * t57 + pkin(7) * t58 + t52;
t18 = -g(2) * t34 + t64;
t43 = cos(qJ(1));
t48 = g(1) * t41 - g(2) * t43;
t30 = t42 * pkin(4) + pkin(3);
t47 = pkin(4) * t59 + t30 * t57 - t34 * t55 + t52;
t46 = (-pkin(3) * t39 - pkin(7) * t38 - pkin(2)) * t64;
t45 = pkin(4) * t56 + t27 + (-t30 * t39 - pkin(2) + t55) * t32;
t36 = qJ(4) + qJ(5);
t35 = t43 * pkin(1);
t33 = cos(t36);
t31 = sin(t36);
t19 = g(1) * t34 + g(2) * t32;
t16 = t18 * t39;
t15 = -g(2) * t58 + t38 * t64;
t14 = t34 * t53 + t59;
t12 = -t32 * t53 + t56;
t10 = t32 * t31 + t33 * t57;
t9 = -t31 * t57 + t32 * t33;
t8 = t34 * t31 - t33 * t60;
t7 = t31 * t60 + t34 * t33;
t6 = -g(1) * t12 - g(2) * t14;
t5 = -g(1) * t11 - g(2) * t13;
t4 = -g(1) * t8 - g(2) * t10;
t3 = -g(1) * t7 - g(2) * t9;
t2 = g(1) * t10 - g(2) * t8 + t33 * t62;
t1 = -g(1) * t9 + g(2) * t7 + t31 * t62;
t17 = [0, 0, 0, 0, 0, 0, t48, g(1) * t43 + g(2) * t41, 0, 0, 0, 0, 0, 0, 0, 0, t18, t19, 0, t48 * pkin(1), 0, 0, 0, 0, 0, 0, t16, -t15, -t19, -g(1) * (t50 - t61) - g(2) * (t35 + t52), 0, 0, 0, 0, 0, 0, t6, t5, t15, -g(1) * (t27 - t61) - g(2) * (t35 + t49) - t46, 0, 0, 0, 0, 0, 0, t4, t3, t15, -g(1) * (t45 - t61) - g(2) * (t35 + t47); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, t19, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t15, -t19, -g(1) * t50 - g(2) * t52, 0, 0, 0, 0, 0, 0, t6, t5, t15, -g(1) * t27 - g(2) * t49 - t46, 0, 0, 0, 0, 0, 0, t4, t3, t15, -g(1) * t45 - g(2) * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, g(1) * t14 - g(2) * t12 + t42 * t62, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t66 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t17;
