% Calculate minimal parameter regressor of gravitation load for
% S5PRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% taug_reg [5x20]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRRP5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:49:17
% EndTime: 2019-12-05 16:49:18
% DurationCPUTime: 0.14s
% Computational Cost: add. (111->37), mult. (166->60), div. (0->0), fcn. (169->8), ass. (0->24)
t14 = sin(pkin(8));
t15 = cos(pkin(8));
t22 = g(1) * t15 + g(2) * t14;
t17 = sin(qJ(2));
t19 = cos(qJ(2));
t3 = -g(3) * t19 + t22 * t17;
t28 = g(3) * t17;
t13 = qJ(3) + qJ(4);
t10 = cos(t13);
t18 = cos(qJ(3));
t7 = t18 * pkin(3) + pkin(4) * t10;
t26 = t14 * t19;
t25 = t15 * t19;
t16 = sin(qJ(3));
t24 = t16 * t19;
t23 = t18 * t19;
t9 = sin(t13);
t1 = -g(1) * (t14 * t10 - t9 * t25) - g(2) * (-t15 * t10 - t9 * t26) + t9 * t28;
t12 = -qJ(5) - pkin(7) - pkin(6);
t6 = -t16 * pkin(3) - pkin(4) * t9;
t5 = pkin(2) + t7;
t4 = t22 * t19 + t28;
t2 = -g(1) * (-t10 * t25 - t14 * t9) - g(2) * (-t10 * t26 + t15 * t9) + t10 * t28;
t8 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, t3, t4, 0, 0, 0, 0, 0, t3 * t18, -t3 * t16, 0, 0, 0, 0, 0, t3 * t10, -t3 * t9, -t4, -g(3) * (-t17 * t12 + t19 * t5) + t22 * (t12 * t19 + t17 * t5); 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t14 * t18 - t15 * t24) - g(2) * (-t14 * t24 - t15 * t18) + t16 * t28, -g(1) * (-t14 * t16 - t15 * t23) - g(2) * (-t14 * t23 + t15 * t16) + t18 * t28, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (t14 * t7 + t6 * t25) - g(2) * (-t15 * t7 + t6 * t26) - t6 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3;];
taug_reg = t8;
