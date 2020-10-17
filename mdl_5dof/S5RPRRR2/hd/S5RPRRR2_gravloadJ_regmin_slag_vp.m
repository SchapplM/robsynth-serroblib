% Calculate minimal parameter regressor of gravitation load for
% S5RPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% taug_reg [5x28]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR2_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR2_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:12:13
% EndTime: 2019-12-05 18:12:14
% DurationCPUTime: 0.09s
% Computational Cost: add. (147->20), mult. (106->30), div. (0->0), fcn. (104->10), ass. (0->18)
t15 = pkin(9) + qJ(3);
t14 = qJ(4) + t15;
t18 = sin(qJ(1));
t19 = cos(qJ(1));
t6 = g(1) * t19 + g(2) * t18;
t5 = g(1) * t18 - g(2) * t19;
t13 = cos(t15);
t12 = sin(t15);
t11 = qJ(5) + t14;
t10 = cos(t14);
t9 = sin(t14);
t8 = cos(t11);
t7 = sin(t11);
t4 = g(3) * t9 + t6 * t10;
t3 = -g(3) * t10 + t6 * t9;
t2 = g(3) * t7 + t6 * t8;
t1 = -g(3) * t8 + t6 * t7;
t16 = [0, t5, t6, t5 * cos(pkin(9)), -t5 * sin(pkin(9)), -t6, -g(1) * (-t18 * pkin(1) + t19 * qJ(2)) - g(2) * (t19 * pkin(1) + t18 * qJ(2)), 0, 0, 0, 0, 0, t5 * t13, -t5 * t12, 0, 0, 0, 0, 0, t5 * t10, -t5 * t9, 0, 0, 0, 0, 0, t5 * t8, -t5 * t7; 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t13 + t6 * t12, g(3) * t12 + t6 * t13, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t16;
