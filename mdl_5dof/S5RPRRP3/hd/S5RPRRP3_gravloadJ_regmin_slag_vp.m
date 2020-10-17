% Calculate minimal parameter regressor of gravitation load for
% S5RPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% taug_reg [5x20]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRP3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:47:42
% EndTime: 2020-01-03 11:47:42
% DurationCPUTime: 0.09s
% Computational Cost: add. (99->25), mult. (89->36), div. (0->0), fcn. (79->8), ass. (0->19)
t16 = cos(qJ(3));
t13 = qJ(3) + qJ(4);
t9 = cos(t13);
t21 = t16 * pkin(3) + pkin(4) * t9;
t12 = qJ(1) + pkin(8);
t6 = sin(t12);
t7 = cos(t12);
t20 = g(2) * t7 + g(3) * t6;
t19 = g(2) * t6 - g(3) * t7;
t15 = sin(qJ(1));
t17 = cos(qJ(1));
t18 = -g(2) * t17 - g(3) * t15;
t8 = sin(t13);
t2 = -g(1) * t9 + t19 * t8;
t14 = sin(qJ(3));
t11 = -qJ(5) - pkin(7) - pkin(6);
t3 = pkin(2) + t21;
t1 = g(1) * t8 + t19 * t9;
t4 = [0, t18, g(2) * t15 - g(3) * t17, t18 * pkin(1), 0, 0, 0, 0, 0, -t20 * t16, t20 * t14, 0, 0, 0, 0, 0, -t20 * t9, t20 * t8, -t19, -g(2) * (t17 * pkin(1) - t6 * t11 + t7 * t3) - g(3) * (t15 * pkin(1) + t7 * t11 + t6 * t3); 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t16 + t19 * t14, g(1) * t14 + t19 * t16, 0, 0, 0, 0, 0, t2, t1, 0, -g(1) * t21 - t19 * (-t14 * pkin(3) - pkin(4) * t8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20;];
taug_reg = t4;
