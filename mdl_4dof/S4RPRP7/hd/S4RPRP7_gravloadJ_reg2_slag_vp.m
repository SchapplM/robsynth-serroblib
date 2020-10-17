% Calculate inertial parameters regressor of gravitation load for
% S4RPRP7
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
% Datum: 2019-12-31 16:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RPRP7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP7_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP7_gravloadJ_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:47:17
% EndTime: 2019-12-31 16:47:17
% DurationCPUTime: 0.11s
% Computational Cost: add. (45->28), mult. (102->31), div. (0->0), fcn. (93->4), ass. (0->16)
t14 = sin(qJ(3));
t16 = cos(qJ(3));
t18 = t14 * pkin(3) - t16 * qJ(4);
t15 = sin(qJ(1));
t17 = cos(qJ(1));
t27 = -g(1) * t15 + g(2) * t17;
t26 = -pkin(1) - pkin(5);
t22 = t17 * pkin(1) + t15 * qJ(2);
t20 = g(2) * (t17 * pkin(5) + t22);
t6 = g(1) * t17 + g(2) * t15;
t9 = t17 * qJ(2);
t4 = t6 * t16;
t3 = t6 * t14;
t2 = -g(3) * t14 - t27 * t16;
t1 = g(3) * t16 - t14 * t27;
t5 = [0, 0, 0, 0, 0, 0, -t27, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t6, -g(1) * (-t15 * pkin(1) + t9) - g(2) * t22, 0, 0, 0, 0, 0, 0, -t3, -t4, -t27, -g(1) * (t26 * t15 + t9) - t20, 0, 0, 0, 0, 0, 0, -t3, -t27, t4, -g(1) * (t18 * t17 + t9) - t20 + (-g(1) * t26 - g(2) * t18) * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, -t1, g(3) * t18 + t27 * (pkin(3) * t16 + qJ(4) * t14); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2;];
taug_reg = t5;
