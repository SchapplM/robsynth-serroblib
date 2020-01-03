% Calculate minimal parameter regressor of gravitation load for
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
% taug_reg [5x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPPR9_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR9_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR9_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPPR9_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t23 = sin(qJ(1));
t26 = cos(qJ(1));
t10 = g(1) * t26 + g(2) * t23;
t22 = sin(qJ(2));
t44 = t10 * t22;
t15 = t22 * qJ(3);
t25 = cos(qJ(2));
t30 = t25 * pkin(2) + t15;
t2 = g(3) * t22 + t10 * t25;
t42 = pkin(2) * t22;
t41 = g(1) * t23;
t37 = g(3) * t25;
t36 = t25 * pkin(3);
t21 = sin(qJ(5));
t35 = t23 * t21;
t24 = cos(qJ(5));
t34 = t23 * t24;
t33 = t25 * t26;
t32 = t26 * t21;
t31 = t26 * t24;
t29 = qJ(3) * t25;
t28 = pkin(2) * t33 + t23 * pkin(6) + (pkin(1) + t15) * t26;
t9 = -g(2) * t26 + t41;
t27 = -pkin(1) - t30;
t18 = t26 * pkin(6);
t13 = t26 * t29;
t11 = t23 * t29;
t8 = t9 * t25;
t7 = t9 * t22;
t6 = t22 * t31 - t35;
t5 = -t22 * t32 - t34;
t4 = -t22 * t34 - t32;
t3 = t22 * t35 - t31;
t1 = -t37 + t44;
t12 = [0, t9, t10, 0, 0, 0, 0, 0, t8, -t7, t8, -t10, t7, -g(1) * t18 - g(2) * t28 - t27 * t41, t7, -t8, t10, -g(1) * (-t26 * qJ(4) + t18) - g(2) * (pkin(3) * t33 + t28) + (-g(1) * (t27 - t36) + g(2) * qJ(4)) * t23, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5; 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, 0, -t2, -g(1) * (-t26 * t42 + t13) - g(2) * (-t23 * t42 + t11) - g(3) * t30, -t2, -t1, 0, -g(1) * t13 - g(2) * t11 - g(3) * (t30 + t36) + (pkin(2) + pkin(3)) * t44, 0, 0, 0, 0, 0, -t2 * t24, t2 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 + g(2) * t3 - t21 * t37, g(1) * t6 - g(2) * t4 - t24 * t37;];
taug_reg = t12;
