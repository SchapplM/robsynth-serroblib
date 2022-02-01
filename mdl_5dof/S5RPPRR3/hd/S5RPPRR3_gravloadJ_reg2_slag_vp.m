% Calculate inertial parameters regressor of gravitation load for
% S5RPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:15
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:14:36
% EndTime: 2022-01-23 09:14:36
% DurationCPUTime: 0.14s
% Computational Cost: add. (163->43), mult. (116->46), div. (0->0), fcn. (105->10), ass. (0->25)
t23 = sin(qJ(1));
t27 = t23 * pkin(1);
t21 = cos(pkin(9));
t9 = t21 * pkin(3) + pkin(2);
t22 = -pkin(6) - qJ(3);
t18 = pkin(9) + qJ(4);
t19 = qJ(1) + pkin(8);
t11 = sin(t19);
t13 = cos(t19);
t4 = g(1) * t13 + g(2) * t11;
t3 = g(1) * t11 - g(2) * t13;
t24 = cos(qJ(1));
t26 = g(1) * t23 - g(2) * t24;
t10 = sin(t18);
t12 = cos(t18);
t25 = -g(3) * t12 + t4 * t10;
t17 = pkin(7) - t22;
t16 = t24 * pkin(1);
t14 = qJ(5) + t18;
t8 = cos(t14);
t7 = sin(t14);
t5 = pkin(4) * t12 + t9;
t2 = g(3) * t7 + t4 * t8;
t1 = -g(3) * t8 + t4 * t7;
t6 = [0, 0, 0, 0, 0, 0, t26, g(1) * t24 + g(2) * t23, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t26 * pkin(1), 0, 0, 0, 0, 0, 0, t3 * t21, -t3 * sin(pkin(9)), -t4, -g(1) * (-t11 * pkin(2) + t13 * qJ(3) - t27) - g(2) * (t13 * pkin(2) + t11 * qJ(3) + t16), 0, 0, 0, 0, 0, 0, t3 * t12, -t3 * t10, -t4, -g(1) * (-t11 * t9 - t13 * t22 - t27) - g(2) * (-t11 * t22 + t13 * t9 + t16), 0, 0, 0, 0, 0, 0, t3 * t8, -t3 * t7, -t4, -g(1) * (-t11 * t5 + t17 * t13 - t27) - g(2) * (t11 * t17 + t13 * t5 + t16); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, g(3) * t10 + t4 * t12, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t25 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t6;
