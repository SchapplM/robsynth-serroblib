% Calculate minimal parameter regressor of gravitation load for
% S6RRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% taug_reg [6x28]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP1_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t31 = cos(qJ(2));
t21 = t31 * pkin(2);
t24 = qJ(2) + pkin(10);
t20 = qJ(4) + t24;
t16 = sin(t20);
t17 = cos(t20);
t30 = cos(qJ(5));
t18 = t30 * pkin(5) + pkin(4);
t25 = -qJ(6) - pkin(9);
t37 = -t16 * t25 + t17 * t18;
t53 = t37 + pkin(3) * cos(t24) + t21;
t29 = sin(qJ(1));
t32 = cos(qJ(1));
t14 = g(1) * t32 + g(2) * t29;
t27 = sin(qJ(5));
t47 = g(3) * t16;
t41 = t32 * t30;
t44 = t29 * t27;
t6 = t17 * t44 + t41;
t42 = t32 * t27;
t43 = t29 * t30;
t8 = -t17 * t42 + t43;
t52 = -g(1) * t8 + g(2) * t6 + t27 * t47;
t3 = -g(3) * t17 + t14 * t16;
t26 = -qJ(3) - pkin(7);
t38 = pkin(5) * t27 + pkin(8) - t26;
t13 = g(1) * t29 - g(2) * t32;
t36 = t16 * t18 + t17 * t25;
t35 = -pkin(1) - t53;
t28 = sin(qJ(2));
t33 = -g(3) * t31 + t14 * t28;
t19 = t21 + pkin(1);
t9 = t17 * t41 + t44;
t7 = -t17 * t43 + t42;
t5 = t13 * t16;
t4 = t14 * t17 + t47;
t2 = t3 * t30;
t1 = t3 * t27;
t10 = [0, t13, t14, 0, 0, 0, 0, 0, t13 * t31, -t13 * t28, -t14, -g(1) * (-t29 * t19 - t32 * t26) - g(2) * (t32 * t19 - t29 * t26) 0, 0, 0, 0, 0, t13 * t17, -t5, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t9, -g(1) * t6 - g(2) * t8, t5 (-g(1) * t38 + g(2) * t35) * t32 + (-g(1) * t35 - g(2) * t38) * t29; 0, 0, 0, 0, 0, 0, 0, 0, t33, g(3) * t28 + t14 * t31, 0, t33 * pkin(2), 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(3) * t53 + t14 * (pkin(3) * sin(t24) + t28 * pkin(2) + t36); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(3) * t37 + t14 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, g(1) * t9 - g(2) * t7 + t30 * t47, 0, t52 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3;];
taug_reg  = t10;
