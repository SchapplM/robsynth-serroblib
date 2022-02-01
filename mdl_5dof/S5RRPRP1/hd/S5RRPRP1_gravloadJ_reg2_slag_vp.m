% Calculate inertial parameters regressor of gravitation load for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:20
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRP1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:19:59
% EndTime: 2022-01-20 10:20:00
% DurationCPUTime: 0.15s
% Computational Cost: add. (213->43), mult. (144->44), div. (0->0), fcn. (127->8), ass. (0->30)
t21 = qJ(1) + qJ(2);
t18 = sin(t21);
t33 = pkin(2) * t18;
t24 = sin(qJ(1));
t32 = t24 * pkin(1);
t17 = pkin(8) + t21;
t13 = sin(t17);
t14 = cos(t17);
t19 = cos(t21);
t15 = pkin(2) * t19;
t31 = t14 * pkin(3) + t13 * pkin(7) + t15;
t25 = cos(qJ(4));
t16 = t25 * pkin(4) + pkin(3);
t22 = -qJ(5) - pkin(7);
t30 = -t13 * t22 + t14 * t16 + t15;
t6 = g(1) * t14 + g(2) * t13;
t5 = g(1) * t13 - g(2) * t14;
t7 = g(1) * t18 - g(2) * t19;
t26 = cos(qJ(1));
t29 = g(1) * t24 - g(2) * t26;
t28 = -t13 * pkin(3) + t14 * pkin(7) - t33;
t27 = -t13 * t16 - t14 * t22 - t33;
t23 = sin(qJ(4));
t1 = -g(3) * t25 + t6 * t23;
t20 = t26 * pkin(1);
t8 = g(1) * t19 + g(2) * t18;
t4 = t5 * t25;
t3 = t5 * t23;
t2 = g(3) * t23 + t6 * t25;
t9 = [0, 0, 0, 0, 0, 0, t29, g(1) * t26 + g(2) * t24, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t29 * pkin(1), 0, 0, 0, 0, 0, 0, t5, t6, 0, -g(1) * (-t32 - t33) - g(2) * (t15 + t20), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (t28 - t32) - g(2) * (t20 + t31), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (t27 - t32) - g(2) * (t20 + t30); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t7 * pkin(2), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * t28 - g(2) * t31, 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * t27 - g(2) * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5;];
taug_reg = t9;
