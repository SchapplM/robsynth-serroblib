% Calculate inertial parameters regressor of gravitation load for
% S4RRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRRR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR3_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR3_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:24:37
% EndTime: 2019-12-31 17:24:37
% DurationCPUTime: 0.12s
% Computational Cost: add. (133->32), mult. (145->48), div. (0->0), fcn. (130->8), ass. (0->24)
t22 = -pkin(6) - pkin(5);
t17 = qJ(2) + qJ(3);
t13 = cos(t17);
t20 = cos(qJ(2));
t15 = t20 * pkin(2);
t25 = pkin(3) * t13 + t15;
t19 = sin(qJ(1));
t21 = cos(qJ(1));
t7 = g(1) * t21 + g(2) * t19;
t24 = g(1) * t19 - g(2) * t21;
t12 = sin(t17);
t3 = -g(3) * t13 + t7 * t12;
t18 = sin(qJ(2));
t23 = -g(3) * t20 + t7 * t18;
t16 = -pkin(7) + t22;
t14 = qJ(4) + t17;
t11 = t15 + pkin(1);
t10 = cos(t14);
t9 = sin(t14);
t5 = pkin(1) + t25;
t4 = g(3) * t12 + t7 * t13;
t2 = g(3) * t9 + t7 * t10;
t1 = -g(3) * t10 + t7 * t9;
t6 = [0, 0, 0, 0, 0, 0, t24, t7, 0, 0, 0, 0, 0, 0, 0, 0, t24 * t20, -t24 * t18, -t7, -g(1) * (-t19 * pkin(1) + t21 * pkin(5)) - g(2) * (t21 * pkin(1) + t19 * pkin(5)), 0, 0, 0, 0, 0, 0, t24 * t13, -t24 * t12, -t7, -g(1) * (-t19 * t11 - t21 * t22) - g(2) * (t21 * t11 - t19 * t22), 0, 0, 0, 0, 0, 0, t24 * t10, -t24 * t9, -t7, -g(1) * (-t21 * t16 - t19 * t5) - g(2) * (-t19 * t16 + t21 * t5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, g(3) * t18 + t7 * t20, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t23 * pkin(2), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(3) * t25 - t7 * (-t18 * pkin(2) - pkin(3) * t12); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t3 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t6;
