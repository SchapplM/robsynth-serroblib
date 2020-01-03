% Calculate minimal parameter regressor of gravitation load for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% taug_reg [5x16]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t19 = qJ(1) + qJ(2);
t14 = pkin(8) + t19;
t10 = cos(t14);
t15 = sin(t19);
t11 = pkin(2) * t15;
t23 = cos(qJ(4));
t13 = t23 * pkin(4) + pkin(3);
t20 = -qJ(5) - pkin(7);
t9 = sin(t14);
t29 = t10 * t20 + t9 * t13 + t11;
t16 = cos(t19);
t12 = pkin(2) * t16;
t28 = t10 * t13 - t9 * t20 + t12;
t27 = g(2) * t9 - g(3) * t10;
t26 = g(2) * t10 + g(3) * t9;
t5 = -g(2) * t16 - g(3) * t15;
t21 = sin(qJ(4));
t25 = -g(1) * t23 + t27 * t21;
t24 = cos(qJ(1));
t22 = sin(qJ(1));
t18 = t24 * pkin(1);
t17 = t22 * pkin(1);
t4 = g(2) * t15 - g(3) * t16;
t2 = t26 * t23;
t1 = t26 * t21;
t3 = [0, -g(2) * t24 - g(3) * t22, g(2) * t22 - g(3) * t24, 0, t5, t4, -g(2) * (t12 + t18) - g(3) * (t11 + t17), 0, 0, 0, 0, 0, -t2, t1, -t27, -g(2) * (t18 + t28) - g(3) * (t17 + t29); 0, 0, 0, 0, t5, t4, t5 * pkin(2), 0, 0, 0, 0, 0, -t2, t1, -t27, -g(2) * t28 - g(3) * t29; 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, -g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, g(1) * t21 + t27 * t23, 0, t25 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26;];
taug_reg = t3;
