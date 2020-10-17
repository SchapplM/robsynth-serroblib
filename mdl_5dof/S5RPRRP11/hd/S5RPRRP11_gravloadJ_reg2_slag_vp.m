% Calculate inertial parameters regressor of gravitation load for
% S5RPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRP11_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP11_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP11_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP11_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:54:32
% EndTime: 2019-12-31 18:54:33
% DurationCPUTime: 0.24s
% Computational Cost: add. (218->65), mult. (301->86), div. (0->0), fcn. (310->8), ass. (0->42)
t31 = sin(qJ(1));
t33 = cos(qJ(1));
t13 = g(1) * t33 + g(2) * t31;
t26 = pkin(8) + qJ(3);
t23 = sin(t26);
t52 = t13 * t23;
t24 = cos(t26);
t41 = t24 * pkin(3) + t23 * pkin(7);
t49 = g(3) * t23;
t48 = t23 * t33;
t47 = t24 * t33;
t30 = sin(qJ(4));
t46 = t31 * t30;
t32 = cos(qJ(4));
t45 = t31 * t32;
t29 = -pkin(6) - qJ(2);
t44 = t33 * t29;
t43 = t33 * t30;
t42 = t33 * t32;
t28 = cos(pkin(8));
t20 = t28 * pkin(2) + pkin(1);
t11 = t33 * t20;
t40 = pkin(3) * t47 + pkin(7) * t48 + t11;
t7 = t24 * t46 + t42;
t9 = t24 * t43 - t45;
t39 = g(1) * t7 - g(2) * t9;
t12 = g(1) * t31 - g(2) * t33;
t38 = pkin(4) * t32 + qJ(5) * t30;
t1 = g(1) * t9 + g(2) * t7 + t30 * t49;
t10 = t24 * t42 + t46;
t8 = t24 * t45 - t43;
t36 = g(1) * t10 + g(2) * t8 + t32 * t49;
t35 = -g(3) * t24 + t52;
t34 = (-g(1) * (-t20 - t41) + g(2) * t29) * t31;
t16 = pkin(7) * t47;
t14 = t31 * t24 * pkin(7);
t6 = t12 * t23;
t5 = t13 * t24 + t49;
t4 = t35 * t32;
t3 = t35 * t30;
t2 = g(1) * t8 - g(2) * t10;
t15 = [0, 0, 0, 0, 0, 0, t12, t13, 0, 0, 0, 0, 0, 0, 0, 0, t12 * t28, -t12 * sin(pkin(8)), -t13, -g(1) * (-t31 * pkin(1) + t33 * qJ(2)) - g(2) * (t33 * pkin(1) + t31 * qJ(2)), 0, 0, 0, 0, 0, 0, t12 * t24, -t6, -t13, -g(1) * (-t31 * t20 - t44) - g(2) * (-t31 * t29 + t11), 0, 0, 0, 0, 0, 0, t2, -t39, t6, g(1) * t44 - g(2) * t40 + t34, 0, 0, 0, 0, 0, 0, t2, t6, t39, -g(1) * (-t8 * pkin(4) - t7 * qJ(5) - t44) - g(2) * (t10 * pkin(4) + t9 * qJ(5) + t40) + t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, t5, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, -t5, -g(1) * (-pkin(3) * t48 + t16) - g(2) * (-t31 * t23 * pkin(3) + t14) - g(3) * t41, 0, 0, 0, 0, 0, 0, t4, -t5, t3, -g(1) * t16 - g(2) * t14 - g(3) * (t38 * t24 + t41) + (pkin(3) + t38) * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t36, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t36, -g(1) * (-t9 * pkin(4) + t10 * qJ(5)) - g(2) * (-t7 * pkin(4) + t8 * qJ(5)) - (-pkin(4) * t30 + qJ(5) * t32) * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t15;
