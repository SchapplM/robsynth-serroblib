% Calculate minimal parameter regressor of gravitation load for
% S5RPRRP2
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
% taug_reg [5x18]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:36
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRP2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:36:27
% EndTime: 2021-01-15 12:36:28
% DurationCPUTime: 0.09s
% Computational Cost: add. (146->28), mult. (102->29), div. (0->0), fcn. (93->8), ass. (0->20)
t15 = qJ(1) + pkin(8);
t14 = qJ(3) + t15;
t11 = sin(t14);
t12 = cos(t14);
t19 = cos(qJ(4));
t13 = t19 * pkin(4) + pkin(3);
t16 = -qJ(5) - pkin(7);
t23 = t11 * t13 + t12 * t16;
t22 = -t11 * t16 + t12 * t13;
t6 = g(2) * t12 + g(3) * t11;
t5 = g(2) * t11 - g(3) * t12;
t18 = sin(qJ(1));
t20 = cos(qJ(1));
t21 = -g(2) * t20 - g(3) * t18;
t17 = sin(qJ(4));
t2 = -g(1) * t19 + t5 * t17;
t4 = t6 * t19;
t3 = t6 * t17;
t1 = g(1) * t17 + t5 * t19;
t7 = [0, t21, g(2) * t18 - g(3) * t20, t21 * pkin(1), 0, -t6, t5, 0, 0, 0, 0, 0, -t4, t3, -t4, t3, -t5, -g(2) * (pkin(2) * cos(t15) + t20 * pkin(1) + t22) - g(3) * (pkin(2) * sin(t15) + t18 * pkin(1) + t23); 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1); 0, 0, 0, 0, 0, -t6, t5, 0, 0, 0, 0, 0, -t4, t3, -t4, t3, -t5, -g(2) * t22 - g(3) * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, t2, t1, 0, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6;];
taug_reg = t7;
