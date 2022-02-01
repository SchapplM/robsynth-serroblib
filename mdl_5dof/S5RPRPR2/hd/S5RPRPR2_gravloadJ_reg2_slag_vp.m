% Calculate inertial parameters regressor of gravitation load for
% S5RPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:19
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:19:07
% EndTime: 2022-01-23 09:19:07
% DurationCPUTime: 0.14s
% Computational Cost: add. (198->42), mult. (116->44), div. (0->0), fcn. (104->10), ass. (0->28)
t23 = qJ(1) + pkin(8);
t20 = qJ(3) + t23;
t13 = sin(t20);
t14 = cos(t20);
t35 = t14 * pkin(3) + t13 * qJ(4);
t19 = cos(t23);
t28 = cos(qJ(1));
t34 = t28 * pkin(1) + pkin(2) * t19;
t33 = -t13 * pkin(3) + t14 * qJ(4);
t25 = cos(pkin(9));
t15 = t25 * pkin(4) + pkin(3);
t26 = -pkin(7) - qJ(4);
t32 = -t13 * t26 + t14 * t15;
t17 = sin(t23);
t27 = sin(qJ(1));
t31 = -t27 * pkin(1) - pkin(2) * t17;
t6 = g(1) * t14 + g(2) * t13;
t5 = g(1) * t13 - g(2) * t14;
t30 = g(1) * t27 - g(2) * t28;
t29 = -t13 * t15 - t14 * t26;
t22 = pkin(9) + qJ(5);
t18 = cos(t22);
t16 = sin(t22);
t4 = t5 * t25;
t3 = t5 * sin(pkin(9));
t2 = t5 * t18;
t1 = t5 * t16;
t7 = [0, 0, 0, 0, 0, 0, t30, g(1) * t28 + g(2) * t27, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t17 - g(2) * t19, g(1) * t19 + g(2) * t17, 0, t30 * pkin(1), 0, 0, 0, 0, 0, 0, t5, t6, 0, -g(1) * t31 - g(2) * t34, 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (t31 + t33) - g(2) * (t34 + t35), 0, 0, 0, 0, 0, 0, t2, -t1, -t6, -g(1) * (t29 + t31) - g(2) * (t32 + t34); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * t33 - g(2) * t35, 0, 0, 0, 0, 0, 0, t2, -t1, -t6, -g(1) * t29 - g(2) * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t18 + t6 * t16, g(3) * t16 + t6 * t18, 0, 0;];
taug_reg = t7;
