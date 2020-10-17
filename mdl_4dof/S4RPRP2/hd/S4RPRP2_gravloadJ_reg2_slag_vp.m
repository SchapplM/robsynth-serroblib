% Calculate inertial parameters regressor of gravitation load for
% S4RPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RPRP2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP2_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP2_gravloadJ_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:15:22
% EndTime: 2019-05-04 19:15:23
% DurationCPUTime: 0.12s
% Computational Cost: add. (51->28), mult. (100->30), div. (0->0), fcn. (108->4), ass. (0->16)
t17 = sin(qJ(3));
t18 = sin(qJ(1));
t23 = t18 * t17;
t20 = cos(qJ(1));
t22 = t20 * t17;
t21 = t20 * pkin(1) + t18 * qJ(2);
t19 = cos(qJ(3));
t5 = -t20 * t19 - t23;
t6 = -t18 * t19 + t22;
t3 = g(1) * t6 - g(2) * t5;
t4 = g(1) * t5 + g(2) * t6;
t14 = t20 * qJ(2);
t12 = t19 * pkin(3) + pkin(2);
t8 = g(1) * t20 + g(2) * t18;
t7 = g(1) * t18 - g(2) * t20;
t1 = [0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, -t8, -g(1) * (-t18 * pkin(1) + t14) - g(2) * t21, 0, 0, 0, 0, 0, 0, -t3, t4, 0, -g(1) * (t14 + (-pkin(1) - pkin(2)) * t18) - g(2) * (t20 * pkin(2) + t21) 0, 0, 0, 0, 0, 0, -t3, t4, 0, -g(1) * (pkin(3) * t22 + t14 + (-pkin(1) - t12) * t18) - g(2) * (pkin(3) * t23 + t20 * t12 + t21); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4, 0, t3 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3);];
taug_reg  = t1;
