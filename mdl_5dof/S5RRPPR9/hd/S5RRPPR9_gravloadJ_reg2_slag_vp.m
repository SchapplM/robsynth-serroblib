% Calculate inertial parameters regressor of gravitation load for
% S5RRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPPR9_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR9_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR9_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPPR9_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t28 = cos(qJ(1));
t25 = sin(qJ(1));
t49 = g(2) * t25;
t10 = g(1) * t28 + t49;
t24 = sin(qJ(2));
t53 = t10 * t24;
t27 = cos(qJ(2));
t2 = g(3) * t24 + t10 * t27;
t52 = pkin(2) + pkin(3);
t51 = g(1) * t25;
t47 = g(3) * t27;
t19 = t27 * pkin(2);
t18 = t27 * pkin(3);
t46 = g(2) * qJ(4);
t45 = t24 * t28;
t23 = sin(qJ(5));
t44 = t25 * t23;
t26 = cos(qJ(5));
t43 = t25 * t26;
t42 = t27 * t28;
t41 = t28 * t23;
t40 = t28 * t26;
t16 = t24 * qJ(3);
t39 = t19 + t16;
t38 = t28 * pkin(1) + t25 * pkin(6);
t37 = qJ(3) * t27;
t36 = pkin(7) + t52;
t35 = t18 + t39;
t34 = -pkin(1) - t16;
t33 = g(1) * t36;
t32 = pkin(2) * t42 + t28 * t16 + t38;
t31 = pkin(3) * t42 + t32;
t9 = -g(2) * t28 + t51;
t20 = t28 * pkin(6);
t30 = g(1) * (-t28 * qJ(4) + t20);
t29 = t34 - t19;
t13 = t28 * t37;
t11 = t25 * t37;
t8 = t9 * t27;
t7 = t9 * t24;
t6 = t24 * t40 - t44;
t5 = -t24 * t41 - t43;
t4 = -t24 * t43 - t41;
t3 = t24 * t44 - t40;
t1 = -t47 + t53;
t12 = [0, 0, 0, 0, 0, 0, t9, t10, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, -t10, -g(1) * (-t25 * pkin(1) + t20) - g(2) * t38, 0, 0, 0, 0, 0, 0, t8, -t10, t7, -g(1) * t20 - g(2) * t32 - t29 * t51, 0, 0, 0, 0, 0, 0, t7, -t8, t10, -t30 - g(2) * t31 + (-g(1) * (t29 - t18) + t46) * t25, 0, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5, t8, -t30 - g(2) * (pkin(4) * t45 + pkin(7) * t42 + t31) + (-g(1) * (-t24 * pkin(4) + t34) + t46 + t27 * t33) * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t2, -g(1) * (-pkin(2) * t45 + t13) - g(2) * (-t25 * t24 * pkin(2) + t11) - g(3) * t39, 0, 0, 0, 0, 0, 0, -t2, -t1, 0, -g(1) * t13 - g(2) * t11 - g(3) * t35 + t52 * t53, 0, 0, 0, 0, 0, 0, -t2 * t26, t2 * t23, t1, -g(1) * (pkin(4) * t42 + t13) - g(2) * (t25 * t27 * pkin(4) + t11) - g(3) * (t27 * pkin(7) + t35) + (-g(3) * pkin(4) + t28 * t33 + t36 * t49) * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 + g(2) * t3 - t23 * t47, g(1) * t6 - g(2) * t4 - t26 * t47, 0, 0;];
taug_reg = t12;
