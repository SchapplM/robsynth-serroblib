% Calculate inertial parameters regressor of gravitation load for
% S5RPRRP7
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
% Datum: 2019-12-31 18:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRP7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:45:37
% EndTime: 2019-12-31 18:45:38
% DurationCPUTime: 0.24s
% Computational Cost: add. (232->61), mult. (285->84), div. (0->0), fcn. (294->8), ass. (0->41)
t26 = qJ(1) + pkin(8);
t21 = sin(t26);
t22 = cos(t26);
t11 = g(1) * t22 + g(2) * t21;
t28 = sin(qJ(3));
t53 = t11 * t28;
t31 = cos(qJ(3));
t44 = t31 * pkin(3) + t28 * pkin(7);
t52 = g(1) * t21;
t49 = g(3) * t28;
t48 = t22 * t28;
t47 = t22 * t31;
t27 = sin(qJ(4));
t46 = t27 * t31;
t30 = cos(qJ(4));
t45 = t30 * t31;
t32 = cos(qJ(1));
t43 = t32 * pkin(1) + t22 * pkin(2) + t21 * pkin(6);
t29 = sin(qJ(1));
t42 = -t29 * pkin(1) + t22 * pkin(6);
t6 = t21 * t46 + t22 * t30;
t8 = -t21 * t30 + t22 * t46;
t41 = g(1) * t6 - g(2) * t8;
t40 = pkin(3) * t47 + pkin(7) * t48 + t43;
t39 = -g(2) * t22 + t52;
t38 = g(1) * t29 - g(2) * t32;
t37 = pkin(4) * t30 + qJ(5) * t27;
t1 = g(1) * t8 + g(2) * t6 + t27 * t49;
t7 = t21 * t45 - t22 * t27;
t9 = t21 * t27 + t22 * t45;
t35 = g(1) * t9 + g(2) * t7 + t30 * t49;
t34 = (-pkin(2) - t44) * t52;
t33 = -g(3) * t31 + t53;
t14 = pkin(7) * t47;
t12 = t21 * t31 * pkin(7);
t10 = t39 * t28;
t5 = t11 * t31 + t49;
t4 = t33 * t30;
t3 = t33 * t27;
t2 = g(1) * t7 - g(2) * t9;
t13 = [0, 0, 0, 0, 0, 0, t38, g(1) * t32 + g(2) * t29, 0, 0, 0, 0, 0, 0, 0, 0, t39, t11, 0, t38 * pkin(1), 0, 0, 0, 0, 0, 0, t39 * t31, -t10, -t11, -g(1) * (-t21 * pkin(2) + t42) - g(2) * t43, 0, 0, 0, 0, 0, 0, t2, -t41, t10, -g(1) * t42 - g(2) * t40 - t34, 0, 0, 0, 0, 0, 0, t2, t10, t41, -g(1) * (-t7 * pkin(4) - t6 * qJ(5) + t42) - g(2) * (t9 * pkin(4) + t8 * qJ(5) + t40) - t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, t5, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, -t5, -g(1) * (-pkin(3) * t48 + t14) - g(2) * (-t21 * t28 * pkin(3) + t12) - g(3) * t44, 0, 0, 0, 0, 0, 0, t4, -t5, t3, -g(1) * t14 - g(2) * t12 - g(3) * (t37 * t31 + t44) + (pkin(3) + t37) * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t35, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t35, -g(1) * (-t8 * pkin(4) + t9 * qJ(5)) - g(2) * (-t6 * pkin(4) + t7 * qJ(5)) - (-pkin(4) * t27 + qJ(5) * t30) * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t13;
