% Calculate minimal parameter regressor of gravitation load for
% S5RRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% taug_reg [5x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPPR8_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR8_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t34 = sin(qJ(1));
t36 = cos(qJ(1));
t17 = g(1) * t36 + g(2) * t34;
t33 = sin(qJ(2));
t54 = t17 * t33;
t24 = t33 * qJ(3);
t35 = cos(qJ(2));
t46 = t35 * pkin(2) + t24;
t52 = pkin(3) * t35;
t51 = g(1) * t34;
t48 = t33 * t36;
t47 = t35 * t36;
t45 = qJ(3) * t35;
t44 = pkin(2) * t47 + t34 * pkin(6) + (pkin(1) + t24) * t36;
t16 = -g(2) * t36 + t51;
t30 = pkin(8) + qJ(5);
t22 = sin(t30);
t23 = cos(t30);
t43 = t22 * t35 - t23 * t33;
t42 = t22 * t33 + t23 * t35;
t31 = sin(pkin(8));
t32 = cos(pkin(8));
t41 = t31 * t35 - t32 * t33;
t40 = t31 * t33 + t32 * t35;
t39 = -pkin(1) - t46;
t1 = t43 * t34;
t3 = t22 * t47 - t23 * t48;
t38 = g(1) * t3 + g(2) * t1 + g(3) * t42;
t2 = t42 * t34;
t4 = t42 * t36;
t37 = g(1) * t4 + g(2) * t2 - g(3) * t43;
t27 = t36 * pkin(6);
t20 = t36 * t45;
t18 = t34 * t45;
t14 = t16 * t35;
t13 = t16 * t33;
t10 = t40 * t36;
t9 = t41 * t36;
t8 = t40 * t34;
t7 = t41 * t34;
t6 = g(3) * t33 + t17 * t35;
t5 = -g(3) * t35 + t54;
t11 = [0, t16, t17, 0, 0, 0, 0, 0, t14, -t13, t14, -t17, t13, -g(1) * t27 - g(2) * t44 - t39 * t51, g(1) * t8 - g(2) * t10, -g(1) * t7 + g(2) * t9, t17, -g(1) * (-qJ(4) * t36 + t27) - g(2) * (pkin(3) * t47 + t44) + (-g(1) * (t39 - t52) + g(2) * qJ(4)) * t34, 0, 0, 0, 0, 0, g(1) * t2 - g(2) * t4, -g(1) * t1 + g(2) * t3; 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, t5, 0, -t6, -g(1) * (-pkin(2) * t48 + t20) - g(2) * (-pkin(2) * t33 * t34 + t18) - g(3) * t46, -g(1) * t9 - g(2) * t7 - g(3) * t40, -g(1) * t10 - g(2) * t8 + g(3) * t41, 0, -g(1) * t20 - g(2) * t18 - g(3) * (t46 + t52) + (pkin(2) + pkin(3)) * t54, 0, 0, 0, 0, 0, -t38, -t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t37;];
taug_reg = t11;
