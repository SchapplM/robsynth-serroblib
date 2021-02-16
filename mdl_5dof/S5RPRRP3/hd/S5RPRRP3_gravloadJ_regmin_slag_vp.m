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
% taug_reg [5x22]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:46
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
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
% StartTime: 2021-01-15 12:46:05
% EndTime: 2021-01-15 12:46:06
% DurationCPUTime: 0.10s
% Computational Cost: add. (137->26), mult. (117->36), div. (0->0), fcn. (107->8), ass. (0->21)
t15 = qJ(3) + qJ(4);
t11 = cos(t15);
t18 = cos(qJ(3));
t23 = t18 * pkin(3) + pkin(4) * t11;
t14 = qJ(1) + pkin(8);
t8 = sin(t14);
t9 = cos(t14);
t22 = g(2) * t9 + g(3) * t8;
t21 = g(2) * t8 - g(3) * t9;
t17 = sin(qJ(1));
t19 = cos(qJ(1));
t20 = -g(2) * t19 - g(3) * t17;
t10 = sin(t15);
t2 = -g(1) * t11 + t21 * t10;
t16 = sin(qJ(3));
t13 = qJ(5) + pkin(7) + pkin(6);
t5 = pkin(2) + t23;
t4 = t22 * t11;
t3 = t22 * t10;
t1 = g(1) * t10 + t21 * t11;
t6 = [0, t20, g(2) * t17 - g(3) * t19, t20 * pkin(1), 0, 0, 0, 0, 0, -t22 * t18, t22 * t16, 0, 0, 0, 0, 0, -t4, t3, -t4, t3, -t21, -g(2) * (t19 * pkin(1) + t8 * t13 + t9 * t5) - g(3) * (t17 * pkin(1) - t9 * t13 + t8 * t5); 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t18 + t21 * t16, g(1) * t16 + t21 * t18, 0, 0, 0, 0, 0, t2, t1, t2, t1, 0, -g(1) * t23 - t21 * (-t16 * pkin(3) - pkin(4) * t10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, t2, t1, 0, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22;];
taug_reg = t6;
