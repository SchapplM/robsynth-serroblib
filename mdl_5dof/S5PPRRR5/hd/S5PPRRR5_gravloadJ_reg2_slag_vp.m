% Calculate inertial parameters regressor of gravitation load for
% S5PPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PPRRR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR5_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRR5_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:35:45
% EndTime: 2019-12-31 17:35:45
% DurationCPUTime: 0.13s
% Computational Cost: add. (138->30), mult. (152->38), div. (0->0), fcn. (182->8), ass. (0->25)
t17 = cos(pkin(8));
t19 = sin(qJ(3));
t21 = cos(qJ(3));
t28 = sin(pkin(8));
t10 = -t17 * t21 - t28 * t19;
t30 = g(2) * t10;
t29 = t17 * t19;
t27 = qJ(3) + qJ(4);
t23 = sin(t27);
t24 = cos(t27);
t8 = -t17 * t24 - t28 * t23;
t9 = t17 * t23 - t28 * t24;
t26 = t8 * pkin(4) - t9 * pkin(7);
t25 = -t9 * pkin(4) - t8 * pkin(7);
t22 = t28 * t21;
t4 = g(1) * t9 - g(2) * t8;
t5 = g(1) * t8 + g(2) * t9;
t20 = cos(qJ(5));
t18 = sin(qJ(5));
t16 = pkin(3) * t22;
t12 = -g(1) * t28 + g(2) * t17;
t11 = -t22 + t29;
t2 = t4 * t20;
t1 = t4 * t18;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t11 - t30, -g(1) * t10 - g(2) * t11, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t5, 0, -g(1) * t16 + (g(1) * t29 - t30) * pkin(3), 0, 0, 0, 0, 0, 0, t2, -t1, t5, -g(1) * (-pkin(3) * t29 + t16 + t25) - g(2) * (t10 * pkin(3) + t26); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t5, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1, t5, -g(1) * t25 - g(2) * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t20 - t5 * t18, -g(3) * t18 - t5 * t20, 0, 0;];
taug_reg = t3;
