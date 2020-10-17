% Calculate minimal parameter regressor of gravitation load for
% S5RRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% taug_reg [5x28]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPR7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR7_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR7_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:17:16
% EndTime: 2019-12-31 21:17:17
% DurationCPUTime: 0.20s
% Computational Cost: add. (196->49), mult. (232->74), div. (0->0), fcn. (236->10), ass. (0->44)
t28 = cos(qJ(2));
t23 = qJ(2) + qJ(3);
t19 = sin(t23);
t20 = cos(t23);
t37 = t20 * pkin(3) + t19 * qJ(4);
t49 = t28 * pkin(2) + t37;
t27 = sin(qJ(1));
t29 = cos(qJ(1));
t34 = g(1) * t29 + g(2) * t27;
t5 = -g(3) * t20 + t34 * t19;
t48 = pkin(3) * t19;
t47 = g(3) * t19;
t22 = pkin(9) + qJ(5);
t17 = sin(t22);
t45 = t27 * t17;
t18 = cos(t22);
t44 = t27 * t18;
t24 = sin(pkin(9));
t43 = t27 * t24;
t25 = cos(pkin(9));
t42 = t27 * t25;
t41 = t29 * t17;
t40 = t29 * t18;
t39 = t29 * t24;
t38 = t29 * t25;
t36 = qJ(4) * t20;
t26 = sin(qJ(2));
t35 = -pkin(2) * t26 - t48;
t33 = g(1) * t27 - g(2) * t29;
t32 = pkin(1) + t49;
t30 = -pkin(7) - pkin(6);
t13 = t29 * t36;
t12 = t27 * t36;
t11 = t33 * t19;
t10 = t20 * t40 + t45;
t9 = -t20 * t41 + t44;
t8 = -t20 * t44 + t41;
t7 = t20 * t45 + t40;
t6 = t34 * t20 + t47;
t4 = t5 * t25;
t3 = t5 * t24;
t2 = t5 * t18;
t1 = t5 * t17;
t14 = [0, t33, t34, 0, 0, 0, 0, 0, t33 * t28, -t33 * t26, 0, 0, 0, 0, 0, t33 * t20, -t11, -g(1) * (-t20 * t42 + t39) - g(2) * (t20 * t38 + t43), -g(1) * (t20 * t43 + t38) - g(2) * (-t20 * t39 + t42), t11, (g(1) * t30 - g(2) * t32) * t29 + (g(1) * t32 + g(2) * t30) * t27, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t10, -g(1) * t7 - g(2) * t9; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t28 + t34 * t26, g(3) * t26 + t34 * t28, 0, 0, 0, 0, 0, t5, t6, t4, -t3, -t6, -g(1) * (t35 * t29 + t13) - g(2) * (t35 * t27 + t12) - g(3) * t49, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, t4, -t3, -t6, -g(1) * (-t29 * t48 + t13) - g(2) * (-t27 * t48 + t12) - g(3) * t37, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 + g(2) * t7 + t17 * t47, g(1) * t10 - g(2) * t8 + t18 * t47;];
taug_reg = t14;
