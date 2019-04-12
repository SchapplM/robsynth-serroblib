% Calculate minimal parameter regressor of gravitation load for
% S6RRPRRR14V3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% 
% Output:
% taug_reg [6x35]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRR14V3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR14V3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_gravloadJ_regmin_slag_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t26 = sin(qJ(4));
t31 = cos(qJ(4));
t33 = cos(qJ(1));
t40 = t33 * t31;
t28 = sin(qJ(1));
t32 = cos(qJ(2));
t45 = t28 * t32;
t14 = t26 * t45 + t40;
t24 = sin(qJ(6));
t29 = cos(qJ(6));
t41 = t33 * t26;
t15 = t31 * t45 - t41;
t30 = cos(qJ(5));
t25 = sin(qJ(5));
t27 = sin(qJ(2));
t49 = t27 * t25;
t4 = t15 * t30 + t28 * t49;
t61 = -t14 * t29 + t4 * t24;
t60 = t14 * t24 + t4 * t29;
t59 = -g(1) * t33 - g(2) * t28;
t10 = -g(3) * t32 - t27 * t59;
t56 = g(3) * t27;
t11 = -t32 * t59 + t56;
t52 = t24 * t30;
t51 = t26 * t27;
t50 = t26 * t32;
t48 = t27 * t30;
t47 = t27 * t31;
t46 = t27 * t33;
t44 = t29 * t30;
t43 = t32 * t25;
t42 = t32 * t30;
t39 = t24 * t51;
t38 = t29 * t51;
t37 = g(1) * t28 - g(2) * t33;
t3 = -t15 * t25 + t28 * t48;
t13 = t30 * t47 - t43;
t36 = t25 * t47 + t42;
t19 = t37 * t27;
t18 = t28 * t26 + t32 * t40;
t6 = -t18 * t25 + t30 * t46;
t35 = g(1) * t6 + g(2) * t3 - g(3) * t36;
t17 = -t28 * t31 + t32 * t41;
t34 = g(1) * t17 + g(2) * t14 + g(3) * t51;
t20 = t37 * t32;
t16 = t31 * t42 + t49;
t9 = t13 * t33;
t8 = t13 * t28;
t7 = t18 * t30 + t25 * t46;
t2 = t17 * t24 + t7 * t29;
t1 = t17 * t29 - t7 * t24;
t5 = [0, t37, -t59, 0, 0, 0, 0, 0, t20, -t19, t20, t59, t19, qJ(3) * t19, 0, 0, 0, 0, 0, g(1) * t15 - g(2) * t18, -g(1) * t14 + g(2) * t17, 0, 0, 0, 0, 0, g(1) * t4 - g(2) * t7, g(1) * t3 - g(2) * t6, 0, 0, 0, 0, 0, g(1) * t60 - g(2) * t2, -g(1) * t61 - g(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, t10, t11, t10, 0, -t11, -t11 * qJ(3), 0, 0, 0, 0, 0, t10 * t31, -t10 * t26, 0, 0, 0, 0, 0, g(1) * t9 + g(2) * t8 - g(3) * t16, -g(3) * (-t31 * t43 + t48) + t59 * t36, 0, 0, 0, 0, 0, -g(1) * (-t9 * t29 - t33 * t39) - g(2) * (-t28 * t39 - t8 * t29) - g(3) * (t16 * t29 + t24 * t50) -g(1) * (t9 * t24 - t33 * t38) - g(2) * (t8 * t24 - t28 * t38) - g(3) * (-t16 * t24 + t29 * t50); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, g(1) * t18 + g(2) * t15 + g(3) * t47, 0, 0, 0, 0, 0, t34 * t30, -t34 * t25, 0, 0, 0, 0, 0, -g(1) * (-t17 * t44 + t18 * t24) - g(2) * (-t14 * t44 + t15 * t24) - (t24 * t31 - t26 * t44) * t56, -g(1) * (t17 * t52 + t18 * t29) - g(2) * (t14 * t52 + t15 * t29) - (t26 * t52 + t29 * t31) * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, g(1) * t7 + g(2) * t4 + g(3) * t13, 0, 0, 0, 0, 0, -t35 * t29, t35 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 + g(2) * t61 - g(3) * (-t13 * t24 + t38) g(1) * t2 + g(2) * t60 - g(3) * (-t13 * t29 - t39);];
taug_reg  = t5;
