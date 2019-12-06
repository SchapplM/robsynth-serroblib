% Calculate minimal parameter regressor of gravitation load for
% S5PRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% taug_reg [5x17]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRPP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t12 = pkin(7) + qJ(2);
t7 = sin(t12);
t9 = cos(t12);
t3 = g(1) * t9 + g(2) * t7;
t2 = g(1) * t7 - g(2) * t9;
t13 = qJ(3) + pkin(8);
t10 = cos(t13);
t8 = sin(t13);
t19 = t10 * pkin(4) + t8 * qJ(5);
t15 = sin(qJ(3));
t16 = cos(qJ(3));
t17 = -g(3) * t16 + t3 * t15;
t14 = -qJ(4) - pkin(6);
t11 = t16 * pkin(3);
t6 = t11 + pkin(2);
t4 = t9 * t6;
t1 = -g(3) * t10 + t3 * t8;
t5 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, -g(3); 0, 0, t2, t3, 0, 0, 0, 0, 0, t2 * t16, -t2 * t15, -t3, -g(1) * (-t9 * t14 - t7 * t6) - g(2) * (-t7 * t14 + t4), t2 * t10, -t3, t2 * t8, -g(2) * t4 + (g(1) * t14 - g(2) * t19) * t9 + (-g(1) * (-t19 - t6) + g(2) * t14) * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, g(3) * t15 + t3 * t16, 0, t17 * pkin(3), t1, 0, -g(3) * t8 - t3 * t10, -g(3) * (t11 + t19) + t3 * (pkin(3) * t15 + pkin(4) * t8 - qJ(5) * t10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t5;
