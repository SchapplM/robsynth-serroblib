% Calculate minimal parameter regressor of gravitation load for
% S5RRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% taug_reg [5x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:28:27
% EndTime: 2019-12-05 18:28:27
% DurationCPUTime: 0.09s
% Computational Cost: add. (141->20), mult. (108->30), div. (0->0), fcn. (103->8), ass. (0->20)
t13 = qJ(2) + pkin(9) + qJ(4);
t16 = sin(qJ(1));
t18 = cos(qJ(1));
t6 = g(1) * t18 + g(2) * t16;
t5 = g(1) * t16 - g(2) * t18;
t15 = sin(qJ(2));
t17 = cos(qJ(2));
t19 = -g(3) * t17 + t6 * t15;
t14 = -qJ(3) - pkin(6);
t12 = t17 * pkin(2) + pkin(1);
t11 = qJ(5) + t13;
t10 = cos(t13);
t9 = sin(t13);
t8 = cos(t11);
t7 = sin(t11);
t4 = g(3) * t9 + t6 * t10;
t3 = -g(3) * t10 + t6 * t9;
t2 = g(3) * t7 + t6 * t8;
t1 = -g(3) * t8 + t6 * t7;
t20 = [0, t5, t6, 0, 0, 0, 0, 0, t5 * t17, -t5 * t15, -t6, -g(1) * (-t16 * t12 - t18 * t14) - g(2) * (t18 * t12 - t16 * t14), 0, 0, 0, 0, 0, t5 * t10, -t5 * t9, 0, 0, 0, 0, 0, t5 * t8, -t5 * t7; 0, 0, 0, 0, 0, 0, 0, 0, t19, g(3) * t15 + t6 * t17, 0, t19 * pkin(2), 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t20;
