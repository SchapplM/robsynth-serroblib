% Calculate minimal parameter regressor of gravitation load for
% S5RRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% taug_reg [5x23]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRP8_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP8_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t27 = sin(qJ(2));
t19 = t27 * qJ(3);
t30 = cos(qJ(2));
t40 = t30 * pkin(2) + t19;
t28 = sin(qJ(1));
t31 = cos(qJ(1));
t12 = g(1) * t31 + g(2) * t28;
t29 = cos(qJ(4));
t26 = sin(qJ(4));
t45 = t27 * t26;
t9 = t30 * t29 + t45;
t50 = g(3) * t9;
t49 = g(1) * t28;
t44 = t27 * t31;
t18 = t29 * pkin(4) + pkin(3);
t43 = t30 * t18;
t42 = t30 * t26;
t41 = t30 * t31;
t39 = qJ(3) * t30;
t38 = pkin(4) * t45;
t37 = t26 * t41;
t36 = pkin(2) * t41 + t28 * pkin(6) + (pkin(1) + t19) * t31;
t11 = -g(2) * t31 + t49;
t35 = -t27 * t29 + t42;
t34 = -pkin(1) - t40;
t3 = t35 * t28;
t5 = -t29 * t44 + t37;
t33 = g(1) * t5 + g(2) * t3 + t50;
t4 = t9 * t28;
t6 = t9 * t31;
t32 = g(1) * t6 + g(2) * t4 - g(3) * t35;
t25 = -qJ(5) - pkin(7);
t22 = t31 * pkin(6);
t16 = t31 * t39;
t14 = t28 * t39;
t8 = t11 * t30;
t7 = t11 * t27;
t2 = g(3) * t27 + t12 * t30;
t1 = -g(3) * t30 + t12 * t27;
t10 = [0, t11, t12, 0, 0, 0, 0, 0, t8, -t7, t8, -t12, t7, -g(1) * t22 - g(2) * t36 - t34 * t49, 0, 0, 0, 0, 0, g(1) * t4 - g(2) * t6, -g(1) * t3 + g(2) * t5, t12, -g(1) * (t31 * t25 + t22) - g(2) * (t18 * t41 + t31 * t38 + t36) + (-g(1) * (t34 - t38 - t43) - g(2) * t25) * t28; 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, 0, -t2, -g(1) * (-pkin(2) * t44 + t16) - g(2) * (-t28 * t27 * pkin(2) + t14) - g(3) * t40, 0, 0, 0, 0, 0, -t33, -t32, 0, -g(1) * (pkin(4) * t37 + t16) - g(2) * (t28 * pkin(4) * t42 + t14) - g(3) * (t40 + t43) + (-g(3) * pkin(4) * t26 + t12 * (pkin(2) + t18)) * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, t32, 0, (t12 * t35 + t50) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11;];
taug_reg = t10;
