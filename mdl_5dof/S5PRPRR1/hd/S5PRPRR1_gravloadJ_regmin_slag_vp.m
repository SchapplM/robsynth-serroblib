% Calculate minimal parameter regressor of gravitation load for
% S5PRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% taug_reg [5x22]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRPRR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR1_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR1_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:43:04
% EndTime: 2019-12-05 15:43:04
% DurationCPUTime: 0.08s
% Computational Cost: add. (105->19), mult. (68->24), div. (0->0), fcn. (66->8), ass. (0->14)
t12 = pkin(9) + qJ(4);
t13 = pkin(8) + qJ(2);
t10 = cos(t13);
t8 = sin(t13);
t3 = g(1) * t8 - g(2) * t10;
t4 = g(1) * t10 + g(2) * t8;
t11 = qJ(5) + t12;
t9 = cos(t12);
t7 = sin(t12);
t6 = cos(t11);
t5 = sin(t11);
t2 = g(3) * t5 + t4 * t6;
t1 = -g(3) * t6 + t4 * t5;
t14 = [-g(3), 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t3, t4, t3 * cos(pkin(9)), -t3 * sin(pkin(9)), -t4, -g(1) * (-t8 * pkin(2) + t10 * qJ(3)) - g(2) * (t10 * pkin(2) + t8 * qJ(3)), 0, 0, 0, 0, 0, t3 * t9, -t3 * t7, 0, 0, 0, 0, 0, t3 * t6, -t3 * t5; 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t9 + t4 * t7, g(3) * t7 + t4 * t9, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t14;
