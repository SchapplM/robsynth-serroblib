% Calculate inertial parameters regressor of gravitation load for
% S4RRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRPP1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP1_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP1_gravloadJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:19:13
% EndTime: 2019-05-04 19:19:13
% DurationCPUTime: 0.08s
% Computational Cost: add. (102->24), mult. (62->28), div. (0->0), fcn. (52->6), ass. (0->20)
t15 = qJ(1) + qJ(2);
t12 = sin(t15);
t22 = pkin(2) * t12;
t16 = sin(qJ(1));
t21 = t16 * pkin(1);
t13 = cos(t15);
t10 = pkin(2) * t13;
t11 = pkin(6) + t15;
t8 = sin(t11);
t9 = cos(t11);
t20 = t9 * pkin(3) + t8 * qJ(4) + t10;
t3 = g(1) * t12 - g(2) * t13;
t17 = cos(qJ(1));
t19 = g(1) * t16 - g(2) * t17;
t18 = -t8 * pkin(3) + t9 * qJ(4) - t22;
t14 = t17 * pkin(1);
t4 = g(1) * t13 + g(2) * t12;
t2 = g(1) * t9 + g(2) * t8;
t1 = g(1) * t8 - g(2) * t9;
t5 = [0, 0, 0, 0, 0, 0, t19, g(1) * t17 + g(2) * t16, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t19 * pkin(1), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (-t21 - t22) - g(2) * (t10 + t14) 0, 0, 0, 0, 0, 0, t1, 0, -t2, -g(1) * (t18 - t21) - g(2) * (t14 + t20); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t3 * pkin(2), 0, 0, 0, 0, 0, 0, t1, 0, -t2, -g(1) * t18 - g(2) * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t5;
