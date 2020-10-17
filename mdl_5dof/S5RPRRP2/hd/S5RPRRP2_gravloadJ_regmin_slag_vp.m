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
% taug_reg [5x16]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:45:29
% EndTime: 2020-01-03 11:45:29
% DurationCPUTime: 0.13s
% Computational Cost: add. (114->26), mult. (76->29), div. (0->0), fcn. (67->8), ass. (0->19)
t13 = qJ(1) + pkin(8);
t12 = qJ(3) + t13;
t10 = cos(t12);
t17 = cos(qJ(4));
t11 = t17 * pkin(4) + pkin(3);
t14 = -qJ(5) - pkin(7);
t9 = sin(t12);
t22 = t10 * t14 + t9 * t11;
t21 = t10 * t11 - t9 * t14;
t3 = g(2) * t9 - g(3) * t10;
t4 = g(2) * t10 + g(3) * t9;
t16 = sin(qJ(1));
t18 = cos(qJ(1));
t20 = -g(2) * t18 - g(3) * t16;
t15 = sin(qJ(4));
t19 = -g(1) * t17 + t3 * t15;
t2 = t4 * t17;
t1 = t4 * t15;
t5 = [0, t20, g(2) * t16 - g(3) * t18, t20 * pkin(1), 0, -t4, t3, 0, 0, 0, 0, 0, -t2, t1, -t3, -g(2) * (pkin(2) * cos(t13) + t18 * pkin(1) + t21) - g(3) * (pkin(2) * sin(t13) + t16 * pkin(1) + t22); 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1); 0, 0, 0, 0, 0, -t4, t3, 0, 0, 0, 0, 0, -t2, t1, -t3, -g(2) * t21 - g(3) * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, g(1) * t15 + t3 * t17, 0, t19 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4;];
taug_reg = t5;
