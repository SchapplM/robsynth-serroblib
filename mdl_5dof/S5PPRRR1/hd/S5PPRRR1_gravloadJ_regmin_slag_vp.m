% Calculate minimal parameter regressor of gravitation load for
% S5PPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% taug_reg [5x15]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PPRRR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:12:51
% EndTime: 2019-12-05 15:12:51
% DurationCPUTime: 0.10s
% Computational Cost: add. (97->20), mult. (84->28), div. (0->0), fcn. (88->8), ass. (0->21)
t11 = sin(pkin(8));
t12 = cos(pkin(8));
t15 = g(1) * t12 + g(2) * t11;
t10 = pkin(9) + qJ(3);
t9 = qJ(4) + t10;
t5 = sin(t9);
t6 = cos(t9);
t3 = -g(3) * t6 + t15 * t5;
t21 = g(3) * t5;
t13 = sin(qJ(5));
t19 = t11 * t13;
t14 = cos(qJ(5));
t18 = t11 * t14;
t17 = t12 * t13;
t16 = t12 * t14;
t8 = cos(t10);
t7 = sin(t10);
t4 = t15 * t6 + t21;
t2 = t3 * t14;
t1 = t3 * t13;
t20 = [-g(3), -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, -g(1) * t11 + g(2) * t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -g(3) * t8 + t15 * t7, g(3) * t7 + t15 * t8, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t6 * t17 + t18) - g(2) * (-t6 * t19 - t16) + t13 * t21, -g(1) * (-t6 * t16 - t19) - g(2) * (-t6 * t18 + t17) + t14 * t21;];
taug_reg = t20;
