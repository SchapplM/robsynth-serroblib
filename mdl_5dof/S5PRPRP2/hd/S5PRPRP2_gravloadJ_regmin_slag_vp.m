% Calculate minimal parameter regressor of gravitation load for
% S5PRPRP2
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
% taug_reg [5x17]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRPRP2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t16 = sin(pkin(8));
t17 = cos(pkin(8));
t20 = cos(qJ(4));
t32 = -t16 * (-qJ(5) - pkin(6)) + (t20 * pkin(4) + pkin(3)) * t17;
t15 = pkin(7) + qJ(2);
t13 = sin(t15);
t14 = cos(t15);
t19 = sin(qJ(4));
t23 = t17 * t19;
t1 = t13 * t23 + t14 * t20;
t28 = g(3) * t16;
t3 = t13 * t20 - t14 * t23;
t31 = -g(1) * t3 + g(2) * t1 + t19 * t28;
t27 = t14 * pkin(2) + t13 * qJ(3);
t25 = t14 * t19;
t22 = t17 * t20;
t7 = g(1) * t14 + g(2) * t13;
t6 = g(1) * t13 - g(2) * t14;
t9 = t14 * qJ(3);
t5 = t6 * t16;
t4 = t13 * t19 + t14 * t22;
t2 = -t13 * t22 + t25;
t8 = [-g(3), 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, t6, t7, t6 * t17, -t5, -t7, -g(1) * (-t13 * pkin(2) + t9) - g(2) * t27, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t4, -g(1) * t1 - g(2) * t3, t5, -g(1) * (pkin(4) * t25 + t9) - g(2) * (t32 * t14 + t27) + (-g(1) * (-pkin(2) - t32) - g(2) * pkin(4) * t19) * t13; 0, 0, 0, 0, 0, 0, 0, -t6, 0, 0, 0, 0, 0, 0, 0, 0, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, g(1) * t4 - g(2) * t2 + t20 * t28, 0, t31 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t17 - t7 * t16;];
taug_reg = t8;
