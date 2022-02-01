% Calculate minimal parameter regressor of gravitation load for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% taug_reg [5x23]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:18
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:17:25
% EndTime: 2022-01-20 11:17:25
% DurationCPUTime: 0.12s
% Computational Cost: add. (185->35), mult. (164->60), div. (0->0), fcn. (188->10), ass. (0->36)
t39 = g(3) * sin(pkin(9));
t26 = qJ(1) + qJ(2);
t22 = sin(t26);
t28 = cos(pkin(9));
t38 = t22 * t28;
t24 = cos(t26);
t37 = t24 * t28;
t29 = sin(qJ(4));
t36 = t28 * t29;
t31 = cos(qJ(4));
t35 = t28 * t31;
t34 = t24 * pkin(2) + t22 * qJ(3);
t33 = -t22 * pkin(2) + t24 * qJ(3);
t16 = g(1) * t22 - g(2) * t24;
t32 = cos(qJ(1));
t30 = sin(qJ(1));
t25 = qJ(4) + qJ(5);
t23 = cos(t25);
t21 = sin(t25);
t17 = g(1) * t24 + g(2) * t22;
t15 = t16 * t28;
t14 = t22 * t29 + t24 * t35;
t13 = t22 * t31 - t24 * t36;
t12 = -t22 * t35 + t24 * t29;
t11 = t22 * t36 + t24 * t31;
t10 = t22 * t21 + t23 * t37;
t9 = -t21 * t37 + t22 * t23;
t8 = t24 * t21 - t23 * t38;
t7 = t21 * t38 + t24 * t23;
t6 = -g(1) * t12 - g(2) * t14;
t5 = -g(1) * t11 - g(2) * t13;
t4 = -g(1) * t8 - g(2) * t10;
t3 = -g(1) * t7 - g(2) * t9;
t2 = g(1) * t10 - g(2) * t8 + t23 * t39;
t1 = -g(1) * t9 + g(2) * t7 + t21 * t39;
t18 = [0, g(1) * t30 - g(2) * t32, g(1) * t32 + g(2) * t30, 0, t16, t17, t15, -t17, -g(1) * (-t30 * pkin(1) + t33) - g(2) * (t32 * pkin(1) + t34), 0, 0, 0, 0, 0, t6, t5, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, t16, t17, t15, -t17, -g(1) * t33 - g(2) * t34, 0, 0, 0, 0, 0, t6, t5, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t13 + g(2) * t11 + t29 * t39, g(1) * t14 - g(2) * t12 + t31 * t39, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t18;
