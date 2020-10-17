% Calculate inertial parameters regressor of gravitation load for
% S5RPRPR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR15_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR15_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR15_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR15_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:37:24
% EndTime: 2019-12-31 18:37:25
% DurationCPUTime: 0.27s
% Computational Cost: add. (125->59), mult. (229->84), div. (0->0), fcn. (227->8), ass. (0->41)
t26 = cos(qJ(1));
t24 = sin(qJ(1));
t49 = g(1) * t24;
t51 = g(2) * t26 - t49;
t23 = sin(qJ(3));
t25 = cos(qJ(3));
t6 = -g(3) * t23 - t25 * t51;
t50 = -pkin(1) - pkin(6);
t46 = g(3) * t25;
t45 = t23 * t26;
t19 = pkin(8) + qJ(5);
t12 = sin(t19);
t44 = t24 * t12;
t13 = cos(t19);
t43 = t24 * t13;
t20 = sin(pkin(8));
t42 = t24 * t20;
t21 = cos(pkin(8));
t41 = t24 * t21;
t22 = -pkin(7) - qJ(4);
t40 = t25 * t22;
t39 = t26 * t12;
t38 = t26 * t13;
t37 = t26 * t20;
t36 = t26 * t21;
t35 = t26 * pkin(1) + t24 * qJ(2);
t34 = t25 * qJ(4);
t33 = t26 * pkin(6) + t35;
t32 = g(2) * t33;
t9 = g(1) * t26 + g(2) * t24;
t30 = t23 * pkin(3) - t34;
t11 = t21 * pkin(4) + pkin(3);
t28 = t23 * t11 + t40;
t15 = t26 * qJ(2);
t7 = t9 * t25;
t5 = -g(2) * t45 + t23 * t49 + t46;
t4 = t23 * t38 - t44;
t3 = t23 * t39 + t43;
t2 = t23 * t43 + t39;
t1 = -t23 * t44 + t38;
t8 = [0, 0, 0, 0, 0, 0, -t51, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, -t9, -g(1) * (-t24 * pkin(1) + t15) - g(2) * t35, 0, 0, 0, 0, 0, 0, -t9 * t23, -t7, -t51, -g(1) * (t50 * t24 + t15) - t32, 0, 0, 0, 0, 0, 0, -g(1) * (t23 * t36 - t42) - g(2) * (t23 * t41 + t37), -g(1) * (-t23 * t37 - t41) - g(2) * (-t23 * t42 + t36), t7, -g(1) * (pkin(3) * t45 - t26 * t34 + t15) - t32 + (-g(1) * t50 - g(2) * t30) * t24, 0, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t2, g(1) * t3 - g(2) * t1, t7, -g(1) * (t11 * t45 + t26 * t40 + t15) - g(2) * (pkin(4) * t37 + t33) + (-g(1) * (-pkin(4) * t20 + t50) - g(2) * t28) * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5, 0, 0, 0, 0, 0, 0, 0, 0, -t6 * t21, t6 * t20, -t5, g(3) * t30 + t51 * (pkin(3) * t25 + qJ(4) * t23), 0, 0, 0, 0, 0, 0, -t6 * t13, t6 * t12, -t5, g(3) * t28 + t51 * (t11 * t25 - t22 * t23); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 - g(2) * t3 + t12 * t46, g(1) * t2 - g(2) * t4 + t13 * t46, 0, 0;];
taug_reg = t8;
