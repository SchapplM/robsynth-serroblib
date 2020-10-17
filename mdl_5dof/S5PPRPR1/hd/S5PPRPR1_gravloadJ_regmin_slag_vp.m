% Calculate minimal parameter regressor of gravitation load for
% S5PPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% 
% Output:
% taug_reg [5x16]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PPRPR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:01:19
% EndTime: 2019-12-05 15:01:20
% DurationCPUTime: 0.14s
% Computational Cost: add. (86->24), mult. (89->34), div. (0->0), fcn. (90->8), ass. (0->16)
t11 = sin(pkin(7));
t13 = cos(pkin(7));
t15 = g(1) * t13 + g(2) * t11;
t9 = pkin(8) + qJ(3);
t5 = sin(t9);
t7 = cos(t9);
t1 = -g(3) * t7 + t15 * t5;
t22 = g(3) * t5;
t18 = t11 * t7;
t17 = t13 * t7;
t8 = pkin(9) + qJ(5);
t6 = cos(t8);
t4 = sin(t8);
t3 = -g(1) * t11 + g(2) * t13;
t2 = t15 * t7 + t22;
t10 = [-g(3), -g(3), 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, t3, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, t1, t2, t1 * cos(pkin(9)), -t1 * sin(pkin(9)), -t2, -g(3) * (t7 * pkin(3) + t5 * qJ(4)) + t15 * (pkin(3) * t5 - qJ(4) * t7), 0, 0, 0, 0, 0, t1 * t6, -t1 * t4; 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t11 * t6 - t4 * t17) - g(2) * (-t13 * t6 - t4 * t18) + t4 * t22, -g(1) * (-t11 * t4 - t6 * t17) - g(2) * (t13 * t4 - t6 * t18) + t6 * t22;];
taug_reg = t10;
