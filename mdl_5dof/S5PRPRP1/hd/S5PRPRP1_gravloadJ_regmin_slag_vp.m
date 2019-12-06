% Calculate minimal parameter regressor of gravitation load for
% S5PRPRP1
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
% taug_reg [5x19]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRPRP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP1_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP1_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t14 = pkin(7) + qJ(2);
t10 = sin(t14);
t12 = cos(t14);
t6 = g(1) * t12 + g(2) * t10;
t5 = g(1) * t10 - g(2) * t12;
t13 = pkin(8) + qJ(4);
t11 = cos(t13);
t9 = sin(t13);
t19 = t11 * pkin(4) + t9 * qJ(5);
t16 = cos(pkin(8));
t18 = t16 * pkin(3) + pkin(2) + t19;
t17 = -pkin(6) - qJ(3);
t4 = t5 * t11;
t3 = t5 * t9;
t2 = g(3) * t9 + t6 * t11;
t1 = -g(3) * t11 + t6 * t9;
t7 = [-g(3), 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, t5, t6, t5 * t16, -t5 * sin(pkin(8)), -t6, -g(1) * (-t10 * pkin(2) + t12 * qJ(3)) - g(2) * (t12 * pkin(2) + t10 * qJ(3)), 0, 0, 0, 0, 0, t4, -t3, t4, -t6, t3, (g(1) * t17 - g(2) * t18) * t12 + (g(1) * t18 + g(2) * t17) * t10; 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, 0, -t2, -g(3) * t19 + t6 * (pkin(4) * t9 - qJ(5) * t11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t7;
