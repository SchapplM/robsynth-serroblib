% Calculate minimal parameter regressor of gravitation load for
% S4PRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% taug_reg [4x15]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4PRRP6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:43
% EndTime: 2019-12-31 16:30:43
% DurationCPUTime: 0.12s
% Computational Cost: add. (56->27), mult. (149->41), div. (0->0), fcn. (157->6), ass. (0->22)
t11 = sin(pkin(6));
t12 = cos(pkin(6));
t27 = -g(1) * t12 - g(2) * t11;
t14 = sin(qJ(2));
t24 = g(3) * t14;
t13 = sin(qJ(3));
t16 = cos(qJ(2));
t23 = t13 * t16;
t15 = cos(qJ(3));
t22 = t15 * t16;
t20 = pkin(3) * t15 + qJ(4) * t13 + pkin(2);
t5 = t11 * t23 + t12 * t15;
t7 = -t11 * t15 + t12 * t23;
t1 = g(1) * t7 + g(2) * t5 + t13 * t24;
t6 = t11 * t22 - t12 * t13;
t8 = t11 * t13 + t12 * t22;
t19 = g(1) * t8 + g(2) * t6 + t15 * t24;
t18 = -g(3) * t16 - t27 * t14;
t4 = -t27 * t16 + t24;
t3 = t18 * t15;
t2 = t18 * t13;
t9 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, t18, t4, 0, 0, 0, 0, 0, t3, -t2, t3, -t4, t2, -g(3) * (t14 * pkin(5) + t20 * t16) + t27 * (pkin(5) * t16 - t20 * t14); 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t19, t1, 0, -t19, -g(1) * (-t7 * pkin(3) + t8 * qJ(4)) - g(2) * (-t5 * pkin(3) + t6 * qJ(4)) - (-pkin(3) * t13 + qJ(4) * t15) * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t9;
