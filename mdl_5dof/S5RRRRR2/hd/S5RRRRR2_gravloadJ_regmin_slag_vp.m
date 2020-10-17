% Calculate minimal parameter regressor of gravitation load for
% S5RRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% 
% Output:
% taug_reg [5x27]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_gravloadJ_regmin_slag_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:53:43
% EndTime: 2019-12-05 18:53:43
% DurationCPUTime: 0.12s
% Computational Cost: add. (168->26), mult. (158->40), div. (0->0), fcn. (170->10), ass. (0->34)
t22 = qJ(1) + qJ(2);
t18 = sin(t22);
t20 = cos(t22);
t16 = g(1) * t20 + g(2) * t18;
t21 = qJ(3) + qJ(4);
t17 = sin(t21);
t19 = cos(t21);
t5 = -g(3) * t19 + t16 * t17;
t34 = g(3) * t17;
t23 = sin(qJ(5));
t32 = t18 * t23;
t26 = cos(qJ(5));
t31 = t18 * t26;
t30 = t20 * t23;
t29 = t20 * t26;
t15 = g(1) * t18 - g(2) * t20;
t28 = cos(qJ(1));
t27 = cos(qJ(3));
t25 = sin(qJ(1));
t24 = sin(qJ(3));
t14 = t15 * t27;
t13 = t15 * t24;
t12 = t19 * t29 + t32;
t11 = -t19 * t30 + t31;
t10 = -t19 * t31 + t30;
t9 = t19 * t32 + t29;
t8 = t15 * t19;
t7 = t15 * t17;
t6 = t16 * t19 + t34;
t4 = t5 * t26;
t3 = t5 * t23;
t2 = -g(1) * t10 - g(2) * t12;
t1 = -g(1) * t9 - g(2) * t11;
t33 = [0, g(1) * t25 - g(2) * t28, g(1) * t28 + g(2) * t25, 0, t15, t16, 0, 0, 0, 0, 0, t14, -t13, 0, 0, 0, 0, 0, t8, -t7, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, t15, t16, 0, 0, 0, 0, 0, t14, -t13, 0, 0, 0, 0, 0, t8, -t7, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t27 + t16 * t24, g(3) * t24 + t16 * t27, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t11 + g(2) * t9 + t23 * t34, g(1) * t12 - g(2) * t10 + t26 * t34;];
taug_reg = t33;
