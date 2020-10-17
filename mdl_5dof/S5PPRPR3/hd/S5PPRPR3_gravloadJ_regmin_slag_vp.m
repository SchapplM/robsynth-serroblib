% Calculate minimal parameter regressor of gravitation load for
% S5PPRPR3
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
% taug_reg [5x13]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PPRPR3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:05:22
% EndTime: 2019-12-05 15:05:23
% DurationCPUTime: 0.12s
% Computational Cost: add. (61->32), mult. (106->57), div. (0->0), fcn. (121->10), ass. (0->26)
t9 = sin(pkin(8));
t27 = g(3) * t9;
t13 = sin(qJ(5));
t26 = t13 * t9;
t15 = cos(qJ(5));
t25 = t15 * t9;
t10 = sin(pkin(7));
t11 = cos(pkin(8));
t24 = t10 * t11;
t14 = sin(qJ(3));
t23 = t10 * t14;
t16 = cos(qJ(3));
t22 = t10 * t16;
t12 = cos(pkin(7));
t21 = t11 * t12;
t20 = t12 * t14;
t19 = t12 * t16;
t8 = qJ(3) + pkin(9);
t6 = sin(t8);
t7 = cos(t8);
t18 = -g(1) * (t10 * t7 - t6 * t21) - g(2) * (-t12 * t7 - t6 * t24) + t6 * t27;
t17 = -g(1) * (-t11 * t20 + t22) - g(2) * (-t11 * t23 - t19) + t14 * t27;
t5 = -g(1) * t10 + g(2) * t12;
t4 = t10 * t6 + t7 * t21;
t2 = -t12 * t6 + t7 * t24;
t1 = [-g(3), -g(3), 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, t5, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, t17, -g(1) * (-t11 * t19 - t23) - g(2) * (-t11 * t22 + t20) + t16 * t27, t17 * pkin(3), 0, 0, 0, 0, 0, t18 * t15, -t18 * t13; 0, 0, 0, 0, 0, g(3) * t11 + (-g(1) * t12 - g(2) * t10) * t9, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t12 * t25 - t4 * t13) - g(2) * (t10 * t25 - t2 * t13) - g(3) * (-t11 * t15 - t7 * t26), -g(1) * (-t12 * t26 - t4 * t15) - g(2) * (-t10 * t26 - t2 * t15) - g(3) * (t11 * t13 - t7 * t25);];
taug_reg = t1;
