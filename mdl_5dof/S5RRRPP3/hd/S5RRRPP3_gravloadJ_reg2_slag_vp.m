% Calculate inertial parameters regressor of gravitation load for
% S5RRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPP3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:53:41
% EndTime: 2019-12-31 20:53:42
% DurationCPUTime: 0.21s
% Computational Cost: add. (224->61), mult. (238->61), div. (0->0), fcn. (218->6), ass. (0->40)
t26 = sin(qJ(3));
t25 = qJ(1) + qJ(2);
t19 = sin(t25);
t20 = cos(t25);
t6 = g(1) * t20 + g(2) * t19;
t50 = t6 * t26;
t49 = g(1) * t19;
t27 = sin(qJ(1));
t46 = t27 * pkin(1);
t28 = cos(qJ(3));
t22 = t28 * pkin(3);
t45 = t19 * t26;
t44 = t20 * t26;
t43 = t20 * t28;
t42 = pkin(3) + qJ(5);
t41 = t20 * pkin(2) + t19 * pkin(7);
t21 = t26 * qJ(4);
t40 = t22 + t21;
t39 = qJ(4) * t28;
t38 = t28 * qJ(5);
t16 = t20 * pkin(7);
t37 = t16 - t46;
t36 = -t19 * pkin(2) + t16;
t35 = -pkin(2) - t21;
t34 = pkin(3) * t43 + t20 * t21 + t41;
t29 = cos(qJ(1));
t33 = g(1) * t27 - g(2) * t29;
t32 = t19 * pkin(4) + t20 * t38 + t34;
t31 = (t35 - t22) * t49;
t30 = (-t42 * t28 + t35) * t49;
t23 = t29 * pkin(1);
t17 = t20 * pkin(4);
t10 = t20 * t39;
t7 = t19 * t39;
t5 = -g(2) * t20 + t49;
t4 = -g(2) * t43 + t28 * t49;
t3 = g(1) * t45 - g(2) * t44;
t2 = g(3) * t26 + t6 * t28;
t1 = -g(3) * t28 + t50;
t8 = [0, 0, 0, 0, 0, 0, t33, g(1) * t29 + g(2) * t27, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t33 * pkin(1), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (t36 - t46) - g(2) * (t23 + t41), 0, 0, 0, 0, 0, 0, -t6, -t4, t3, -g(1) * t37 - g(2) * (t23 + t34) - t31, 0, 0, 0, 0, 0, 0, -t6, t3, t4, -g(1) * (t17 + t37) - g(2) * (t23 + t32) - t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * t36 - g(2) * t41, 0, 0, 0, 0, 0, 0, -t6, -t4, t3, -g(1) * t16 - g(2) * t34 - t31, 0, 0, 0, 0, 0, 0, -t6, t3, t4, -g(1) * (t16 + t17) - g(2) * t32 - t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, -t2, -g(1) * (-pkin(3) * t44 + t10) - g(2) * (-pkin(3) * t45 + t7) - g(3) * t40, 0, 0, 0, 0, 0, 0, 0, -t2, t1, -g(1) * t10 - g(2) * t7 - g(3) * (t38 + t40) + t42 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2;];
taug_reg = t8;
