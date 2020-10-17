% Calculate inertial parameters regressor of gravitation load for
% S5RRRPP2
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
% Datum: 2019-12-31 20:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPP2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP2_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP2_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:51:58
% EndTime: 2019-12-31 20:51:59
% DurationCPUTime: 0.23s
% Computational Cost: add. (223->53), mult. (235->60), div. (0->0), fcn. (215->6), ass. (0->38)
t23 = sin(qJ(3));
t22 = qJ(1) + qJ(2);
t16 = sin(t22);
t17 = cos(t22);
t6 = g(1) * t17 + g(2) * t16;
t47 = t6 * t23;
t18 = t23 * qJ(4);
t25 = cos(qJ(3));
t36 = t25 * pkin(3) + t18;
t45 = g(1) * t16;
t24 = sin(qJ(1));
t42 = t24 * pkin(1);
t41 = t25 * pkin(4);
t40 = t16 * t23;
t39 = t17 * t23;
t38 = t17 * t25;
t37 = t17 * pkin(2) + t16 * pkin(7);
t35 = qJ(4) * t25;
t14 = t17 * pkin(7);
t34 = -t16 * pkin(2) + t14;
t33 = pkin(3) * t38 + t17 * t18 + t37;
t32 = -t17 * qJ(5) + t14;
t26 = cos(qJ(1));
t20 = t26 * pkin(1);
t31 = t20 + t33;
t5 = -g(2) * t17 + t45;
t30 = g(1) * t24 - g(2) * t26;
t29 = -pkin(2) - t36;
t28 = t29 * t45;
t27 = (-g(1) * (t29 - t41) + g(2) * qJ(5)) * t16;
t10 = pkin(4) * t38;
t9 = t17 * t35;
t7 = t16 * t35;
t4 = t5 * t25;
t3 = g(1) * t40 - g(2) * t39;
t2 = g(3) * t23 + t6 * t25;
t1 = -g(3) * t25 + t47;
t8 = [0, 0, 0, 0, 0, 0, t30, g(1) * t26 + g(2) * t24, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t30 * pkin(1), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (t34 - t42) - g(2) * (t20 + t37), 0, 0, 0, 0, 0, 0, t4, -t6, t3, -g(1) * (t14 - t42) - g(2) * t31 - t28, 0, 0, 0, 0, 0, 0, t4, t3, t6, -g(1) * (t32 - t42) - g(2) * (t10 + t31) + t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * t34 - g(2) * t37, 0, 0, 0, 0, 0, 0, t4, -t6, t3, -g(1) * t14 - g(2) * t33 - t28, 0, 0, 0, 0, 0, 0, t4, t3, t6, -g(1) * t32 - g(2) * (t10 + t33) + t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t2, -g(1) * (-pkin(3) * t39 + t9) - g(2) * (-pkin(3) * t40 + t7) - g(3) * t36, 0, 0, 0, 0, 0, 0, t1, -t2, 0, -g(1) * t9 - g(2) * t7 - g(3) * (t36 + t41) + (pkin(3) + pkin(4)) * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5;];
taug_reg = t8;
