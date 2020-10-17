% Calculate minimal parameter regressor of gravitation load for
% S4RRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% taug_reg [4x17]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRRP3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP3_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP3_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:14:15
% EndTime: 2019-12-31 17:14:15
% DurationCPUTime: 0.13s
% Computational Cost: add. (105->24), mult. (116->30), div. (0->0), fcn. (109->6), ass. (0->20)
t16 = sin(qJ(3));
t18 = cos(qJ(3));
t22 = t18 * pkin(3) + t16 * qJ(4);
t30 = -pkin(2) - t22;
t15 = qJ(1) + qJ(2);
t14 = cos(t15);
t13 = sin(t15);
t29 = g(1) * t13;
t5 = -g(2) * t14 + t29;
t6 = g(1) * t14 + g(2) * t13;
t23 = t13 * pkin(6) - t30 * t14;
t20 = t30 * t29;
t19 = cos(qJ(1));
t17 = sin(qJ(1));
t11 = t14 * pkin(6);
t4 = t5 * t18;
t3 = t5 * t16;
t2 = g(3) * t16 + t6 * t18;
t1 = -g(3) * t18 + t6 * t16;
t7 = [0, g(1) * t17 - g(2) * t19, g(1) * t19 + g(2) * t17, 0, t5, t6, 0, 0, 0, 0, 0, t4, -t3, t4, -t6, t3, -g(1) * (-t17 * pkin(1) + t11) - g(2) * (t19 * pkin(1) + t23) - t20; 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t4, -t3, t4, -t6, t3, -g(1) * t11 - g(2) * t23 - t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, 0, -t2, -g(3) * t22 + t6 * (pkin(3) * t16 - qJ(4) * t18); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t7;
