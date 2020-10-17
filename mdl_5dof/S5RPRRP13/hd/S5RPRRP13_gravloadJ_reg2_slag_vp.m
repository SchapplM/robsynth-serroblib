% Calculate inertial parameters regressor of gravitation load for
% S5RPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRP13_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP13_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP13_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP13_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:59:37
% EndTime: 2019-12-31 18:59:37
% DurationCPUTime: 0.25s
% Computational Cost: add. (120->66), mult. (295->82), div. (0->0), fcn. (304->6), ass. (0->43)
t27 = sin(qJ(1));
t30 = cos(qJ(1));
t51 = g(2) * t30;
t11 = g(1) * t27 - t51;
t29 = cos(qJ(3));
t26 = sin(qJ(3));
t50 = g(3) * t26;
t53 = t11 * t29 - t50;
t52 = pkin(7) * t26;
t49 = g(3) * t29;
t21 = t29 * pkin(7);
t48 = t26 * t27;
t47 = t26 * t30;
t25 = sin(qJ(4));
t46 = t27 * t25;
t28 = cos(qJ(4));
t45 = t27 * t28;
t44 = t27 * t29;
t43 = t30 * t25;
t42 = t30 * t28;
t41 = pkin(3) * t44 + pkin(7) * t48;
t40 = t30 * pkin(1) + t27 * qJ(2);
t39 = t30 * pkin(6) + t40;
t6 = t26 * t46 - t42;
t8 = t26 * t43 + t45;
t38 = g(1) * t8 + g(2) * t6;
t20 = t30 * qJ(2);
t37 = t20 + (-pkin(1) - pkin(6)) * t27;
t12 = g(1) * t30 + g(2) * t27;
t36 = pkin(4) * t28 + qJ(5) * t25;
t35 = -pkin(3) - t36;
t34 = pkin(3) * t48 - pkin(7) * t44 + t39;
t1 = g(1) * t6 - g(2) * t8 + t25 * t49;
t7 = t26 * t45 + t43;
t9 = t26 * t42 - t46;
t33 = g(1) * t7 - g(2) * t9 + t28 * t49;
t32 = pkin(3) * t47 - t30 * t21 + t37;
t10 = t12 * t29;
t5 = g(1) * t48 - g(2) * t47 + t49;
t4 = t53 * t28;
t3 = t53 * t25;
t2 = -g(1) * t9 - g(2) * t7;
t13 = [0, 0, 0, 0, 0, 0, t11, t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t12, -g(1) * (-t27 * pkin(1) + t20) - g(2) * t40, 0, 0, 0, 0, 0, 0, -t12 * t26, -t10, t11, -g(1) * t37 - g(2) * t39, 0, 0, 0, 0, 0, 0, t2, t38, t10, -g(1) * t32 - g(2) * t34, 0, 0, 0, 0, 0, 0, t2, t10, -t38, -g(1) * (t9 * pkin(4) + t8 * qJ(5) + t32) - g(2) * (t7 * pkin(4) + t6 * qJ(5) + t34); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, t5, 0, 0, 0, 0, 0, 0, 0, 0, -t4, t3, -t5, -g(1) * t41 - g(3) * (-t26 * pkin(3) + t21) - (-pkin(3) * t29 - t52) * t51, 0, 0, 0, 0, 0, 0, -t4, -t5, -t3, -g(1) * (t36 * t44 + t41) - g(3) * t21 - t35 * t50 - (t29 * t35 - t52) * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t33, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t33, -g(1) * (-t6 * pkin(4) + t7 * qJ(5)) - g(2) * (t8 * pkin(4) - t9 * qJ(5)) - (-pkin(4) * t25 + qJ(5) * t28) * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t13;
