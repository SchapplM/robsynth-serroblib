% Calculate inertial parameters regressor of gravitation load for
% S5RPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRR9_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR9_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:48
% EndTime: 2019-12-31 18:02:48
% DurationCPUTime: 0.19s
% Computational Cost: add. (127->47), mult. (268->66), div. (0->0), fcn. (319->8), ass. (0->32)
t35 = sin(pkin(8));
t36 = cos(pkin(8));
t40 = sin(qJ(1));
t41 = cos(qJ(1));
t7 = -t40 * t35 - t41 * t36;
t8 = t41 * t35 - t40 * t36;
t32 = g(1) * t7 + g(2) * t8;
t20 = sin(qJ(4));
t22 = cos(qJ(4));
t46 = -g(3) * t22 + t32 * t20;
t43 = g(3) * t20;
t19 = sin(qJ(5));
t39 = t19 * t22;
t21 = cos(qJ(5));
t38 = t21 * t22;
t37 = t41 * pkin(1) + t40 * qJ(2);
t34 = t41 * pkin(2) + t37;
t33 = g(1) * t8 - g(2) * t7;
t31 = -t40 * pkin(1) + t41 * qJ(2);
t30 = t22 * pkin(4) + t20 * pkin(7);
t28 = t7 * t19 + t8 * t38;
t27 = -t7 * t21 + t8 * t39;
t26 = -t7 * pkin(3) + t8 * pkin(6) + t34;
t25 = -t40 * pkin(2) + t31;
t23 = t8 * pkin(3) + t7 * pkin(6) + t25;
t10 = g(1) * t41 + g(2) * t40;
t9 = g(1) * t40 - g(2) * t41;
t4 = t33 * t20;
t3 = t8 * t19 - t7 * t38;
t2 = t8 * t21 + t7 * t39;
t1 = -t32 * t22 - t43;
t5 = [0, 0, 0, 0, 0, 0, t9, t10, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, -t10, -g(1) * t31 - g(2) * t37, 0, 0, 0, 0, 0, 0, -t33, t32, 0, -g(1) * t25 - g(2) * t34, 0, 0, 0, 0, 0, 0, -t33 * t22, t4, -t32, -g(1) * t23 - g(2) * t26, 0, 0, 0, 0, 0, 0, -g(1) * t28 - g(2) * t3, g(1) * t27 - g(2) * t2, -t4, -g(1) * (t30 * t8 + t23) - g(2) * (-t30 * t7 + t26); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, t1, 0, 0, 0, 0, 0, 0, 0, 0, -t46 * t21, t46 * t19, -t1, g(3) * t30 - t32 * (pkin(4) * t20 - pkin(7) * t22); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t27 - t19 * t43, g(1) * t3 - g(2) * t28 - t21 * t43, 0, 0;];
taug_reg = t5;
