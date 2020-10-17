% Calculate minimal parameter regressor of gravitation load for
% S5RPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
% 
% Output:
% taug_reg [5x23]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:38:15
% EndTime: 2019-12-05 17:38:15
% DurationCPUTime: 0.08s
% Computational Cost: add. (49->20), mult. (78->26), div. (0->0), fcn. (74->6), ass. (0->14)
t13 = sin(qJ(1));
t15 = cos(qJ(1));
t16 = t15 * pkin(1) + t13 * qJ(2);
t4 = g(1) * t15 + g(2) * t13;
t3 = g(1) * t13 - g(2) * t15;
t14 = cos(qJ(4));
t12 = sin(qJ(4));
t11 = qJ(4) + qJ(5);
t8 = t15 * qJ(2);
t6 = cos(t11);
t5 = sin(t11);
t2 = g(3) * t5 - t4 * t6;
t1 = g(3) * t6 + t4 * t5;
t7 = [0, t3, t4, -t3, -t4, -g(1) * (-t13 * pkin(1) + t8) - g(2) * t16, -t4, t3, -g(1) * (t8 + (-pkin(1) - qJ(3)) * t13) - g(2) * (t15 * qJ(3) + t16), 0, 0, 0, 0, 0, t3 * t12, t3 * t14, 0, 0, 0, 0, 0, t3 * t5, t3 * t6; 0, 0, 0, 0, 0, -t3, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t12 - t4 * t14, g(3) * t14 + t4 * t12, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1;];
taug_reg = t7;
