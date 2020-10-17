% Calculate minimal parameter regressor of gravitation load for
% S4RRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% taug_reg [4x17]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRPR4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:02:33
% EndTime: 2019-12-31 17:02:33
% DurationCPUTime: 0.09s
% Computational Cost: add. (85->21), mult. (74->26), div. (0->0), fcn. (70->8), ass. (0->17)
t15 = qJ(1) + qJ(2);
t12 = sin(t15);
t13 = cos(t15);
t21 = t13 * pkin(2) + t12 * qJ(3);
t20 = -t12 * pkin(2) + t13 * qJ(3);
t6 = g(1) * t13 + g(2) * t12;
t5 = g(1) * t12 - g(2) * t13;
t19 = cos(qJ(1));
t18 = sin(qJ(1));
t14 = pkin(7) + qJ(4);
t11 = cos(t14);
t10 = sin(t14);
t4 = t5 * cos(pkin(7));
t3 = t5 * sin(pkin(7));
t2 = t5 * t11;
t1 = t5 * t10;
t7 = [0, g(1) * t18 - g(2) * t19, g(1) * t19 + g(2) * t18, 0, t5, t6, t4, -t3, -t6, -g(1) * (-t18 * pkin(1) + t20) - g(2) * (t19 * pkin(1) + t21), 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, t5, t6, t4, -t3, -t6, -g(1) * t20 - g(2) * t21, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t11 + t6 * t10, g(3) * t10 + t6 * t11;];
taug_reg = t7;
