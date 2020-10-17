% Calculate inertial parameters regressor of gravitation load for
% S4RPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RPRR8_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR8_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR8_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR8_gravloadJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:55:15
% EndTime: 2019-12-31 16:55:15
% DurationCPUTime: 0.10s
% Computational Cost: add. (60->30), mult. (98->34), div. (0->0), fcn. (89->6), ass. (0->17)
t13 = sin(qJ(1));
t15 = cos(qJ(1));
t19 = t15 * pkin(1) + t13 * qJ(2);
t12 = sin(qJ(3));
t18 = pkin(3) * t12;
t4 = g(1) * t15 + g(2) * t13;
t3 = g(1) * t13 - g(2) * t15;
t14 = cos(qJ(3));
t17 = g(3) * t12 - t3 * t14;
t16 = -pkin(6) - pkin(5);
t11 = qJ(3) + qJ(4);
t8 = t15 * qJ(2);
t6 = cos(t11);
t5 = sin(t11);
t2 = g(3) * t5 - t3 * t6;
t1 = g(3) * t6 + t3 * t5;
t7 = [0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, -g(1) * (-t13 * pkin(1) + t8) - g(2) * t19, 0, 0, 0, 0, 0, 0, -t4 * t12, -t4 * t14, t3, -g(1) * (t8 + (-pkin(1) - pkin(5)) * t13) - g(2) * (t15 * pkin(5) + t19), 0, 0, 0, 0, 0, 0, -t4 * t5, -t4 * t6, t3, -g(1) * (t15 * t18 + t8 + (-pkin(1) + t16) * t13) - g(2) * (t13 * t18 - t15 * t16 + t19); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, g(3) * t14 + t3 * t12, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, t17 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0;];
taug_reg = t7;
