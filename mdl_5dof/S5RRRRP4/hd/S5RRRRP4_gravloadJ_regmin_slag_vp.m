% Calculate minimal parameter regressor of gravitation load for
% S5RRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% taug_reg [5x24]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRP4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP4_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP4_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:51:13
% EndTime: 2019-12-31 21:51:13
% DurationCPUTime: 0.18s
% Computational Cost: add. (226->40), mult. (184->52), div. (0->0), fcn. (174->8), ass. (0->32)
t27 = cos(qJ(3));
t23 = qJ(3) + qJ(4);
t18 = sin(t23);
t20 = cos(t23);
t34 = t20 * pkin(4) + t18 * qJ(5);
t39 = t27 * pkin(3) + t34;
t38 = pkin(2) + t39;
t24 = qJ(1) + qJ(2);
t19 = sin(t24);
t37 = t18 * t19;
t21 = cos(t24);
t36 = t18 * t21;
t29 = -pkin(8) - pkin(7);
t35 = t21 * t29;
t33 = qJ(5) * t20;
t32 = t38 * t21;
t25 = sin(qJ(3));
t31 = -pkin(3) * t25 - pkin(4) * t18;
t8 = g(1) * t21 + g(2) * t19;
t7 = g(1) * t19 - g(2) * t21;
t30 = (g(1) * t38 + g(2) * t29) * t19;
t28 = cos(qJ(1));
t26 = sin(qJ(1));
t11 = t21 * t33;
t9 = t19 * t33;
t6 = t7 * t27;
t5 = t7 * t25;
t4 = t7 * t20;
t3 = g(1) * t37 - g(2) * t36;
t2 = g(3) * t18 + t8 * t20;
t1 = -g(3) * t20 + t8 * t18;
t10 = [0, g(1) * t26 - g(2) * t28, g(1) * t28 + g(2) * t26, 0, t7, t8, 0, 0, 0, 0, 0, t6, -t5, 0, 0, 0, 0, 0, t4, -t3, t4, -t8, t3, -g(1) * (-t26 * pkin(1) - t35) - g(2) * (t28 * pkin(1) + t32) + t30; 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, t6, -t5, 0, 0, 0, 0, 0, t4, -t3, t4, -t8, t3, g(1) * t35 - g(2) * t32 + t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t27 + t8 * t25, g(3) * t25 + t8 * t27, 0, 0, 0, 0, 0, t1, t2, t1, 0, -t2, -g(1) * (t31 * t21 + t11) - g(2) * (t31 * t19 + t9) - g(3) * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, 0, -t2, -g(1) * (-pkin(4) * t36 + t11) - g(2) * (-pkin(4) * t37 + t9) - g(3) * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t10;
