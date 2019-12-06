% Calculate minimal parameter regressor of gravitation load for
% S5PRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
% 
% Output:
% taug_reg [5x16]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRPPR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR2_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR2_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t11 = cos(pkin(7));
t9 = sin(pkin(7));
t18 = g(1) * t11 + g(2) * t9;
t7 = qJ(2) + pkin(8);
t3 = sin(t7);
t5 = cos(t7);
t16 = -g(3) * t5 + t18 * t3;
t24 = g(3) * t3;
t6 = pkin(9) + qJ(5);
t2 = sin(t6);
t22 = t9 * t2;
t4 = cos(t6);
t21 = t9 * t4;
t19 = t11 * t5;
t12 = sin(qJ(2));
t13 = cos(qJ(2));
t14 = -g(3) * t13 + t18 * t12;
t1 = -g(1) * t9 + g(2) * t11;
t8 = [-g(3), 0, 0, 0, -g(3), 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, t14, g(3) * t12 + t18 * t13, t14 * pkin(2), t16 * cos(pkin(9)), -t16 * sin(pkin(9)), -t18 * t5 - t24, -g(3) * (t13 * pkin(2) + t5 * pkin(3) + t3 * qJ(4)) + t18 * (pkin(2) * t12 + pkin(3) * t3 - qJ(4) * t5), 0, 0, 0, 0, 0, t16 * t4, -t16 * t2; 0, 0, 0, 0, t1, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t2 * t19 + t21) - g(2) * (-t11 * t4 - t5 * t22) + t2 * t24, -g(1) * (-t4 * t19 - t22) - g(2) * (t11 * t2 - t5 * t21) + t4 * t24;];
taug_reg = t8;
