% Calculate minimal parameter regressor of gravitation load for
% S5PRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% taug_reg [5x14]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRPRP3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t6 = sin(pkin(7));
t7 = cos(pkin(7));
t17 = g(1) * t7 + g(2) * t6;
t5 = qJ(2) + pkin(8);
t3 = sin(t5);
t4 = cos(t5);
t26 = -g(3) * t4 + t17 * t3;
t23 = g(3) * t3;
t9 = sin(qJ(4));
t21 = t6 * t9;
t20 = t7 * t9;
t11 = cos(qJ(4));
t19 = t6 * t11;
t18 = t7 * t11;
t10 = sin(qJ(2));
t12 = cos(qJ(2));
t14 = -g(3) * t12 + t17 * t10;
t13 = -g(1) * (-t4 * t20 + t19) - g(2) * (-t4 * t21 - t18) + t9 * t23;
t8 = -qJ(5) - pkin(6);
t2 = t11 * pkin(4) + pkin(3);
t1 = -g(1) * t6 + g(2) * t7;
t15 = [-g(3), 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, t14, g(3) * t10 + t17 * t12, t14 * pkin(2), 0, 0, 0, 0, 0, t26 * t11, -t26 * t9, -t17 * t4 - t23, -g(3) * (t12 * pkin(2) + t4 * t2 - t3 * t8) + t17 * (pkin(2) * t10 + t2 * t3 + t4 * t8); 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -g(1) * (-t4 * t18 - t21) - g(2) * (-t4 * t19 + t20) + t11 * t23, 0, t13 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26;];
taug_reg = t15;
