% Calculate minimal parameter regressor of gravitation load for
% S5RPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% taug_reg [5x24]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRP13_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP13_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP13_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP13_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t20 = sin(qJ(1));
t23 = cos(qJ(1));
t43 = -g(1) * t20 + g(2) * t23;
t19 = sin(qJ(3));
t22 = cos(qJ(3));
t33 = t22 * pkin(7);
t42 = t19 * pkin(3) - t33;
t41 = -g(3) * t19 - t43 * t22;
t35 = g(3) * t22;
t18 = sin(qJ(4));
t32 = t20 * t18;
t21 = cos(qJ(4));
t31 = t20 * t21;
t30 = t23 * t18;
t29 = t23 * t21;
t28 = t23 * pkin(1) + t20 * qJ(2);
t6 = t19 * t32 - t29;
t8 = t19 * t30 + t31;
t27 = g(1) * t8 + g(2) * t6;
t12 = g(1) * t23 + g(2) * t20;
t26 = pkin(4) * t21 + qJ(5) * t18 + pkin(3);
t1 = g(1) * t6 - g(2) * t8 + t18 * t35;
t7 = t19 * t31 + t30;
t9 = t19 * t29 - t32;
t25 = g(1) * t7 - g(2) * t9 + t21 * t35;
t15 = t23 * qJ(2);
t10 = t12 * t22;
t5 = -t19 * t43 + t35;
t4 = t41 * t21;
t3 = t41 * t18;
t2 = -g(1) * t9 - g(2) * t7;
t11 = [0, -t43, t12, t43, -t12, -g(1) * (-t20 * pkin(1) + t15) - g(2) * t28, 0, 0, 0, 0, 0, -t12 * t19, -t10, 0, 0, 0, 0, 0, t2, t27, t2, t10, -t27, -g(1) * (t9 * pkin(4) + t8 * qJ(5) + t42 * t23 + t15) - g(2) * (t7 * pkin(4) + t23 * pkin(6) + t6 * qJ(5) + t28) + (-g(1) * (-pkin(1) - pkin(6)) - g(2) * t42) * t20; 0, 0, 0, 0, 0, t43, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, t5, 0, 0, 0, 0, 0, -t4, t3, -t4, -t5, -t3, -g(3) * (-t26 * t19 + t33) + t43 * (pkin(7) * t19 + t26 * t22); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t25, t1, 0, -t25, -g(1) * (-t6 * pkin(4) + t7 * qJ(5)) - g(2) * (t8 * pkin(4) - t9 * qJ(5)) - (-pkin(4) * t18 + qJ(5) * t21) * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t11;
