% Calculate inertial parameters regressor of gravitation load for
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
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRPR4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:02:33
% EndTime: 2019-12-31 17:02:34
% DurationCPUTime: 0.11s
% Computational Cost: add. (119->32), mult. (102->36), div. (0->0), fcn. (92->8), ass. (0->25)
t23 = sin(qJ(1));
t30 = t23 * pkin(1);
t19 = qJ(1) + qJ(2);
t15 = sin(t19);
t16 = cos(t19);
t29 = t16 * pkin(2) + t15 * qJ(3);
t28 = -t15 * pkin(2) + t16 * qJ(3);
t21 = cos(pkin(7));
t10 = t21 * pkin(3) + pkin(2);
t22 = -pkin(6) - qJ(3);
t27 = t16 * t10 - t15 * t22;
t6 = g(1) * t16 + g(2) * t15;
t5 = g(1) * t15 - g(2) * t16;
t24 = cos(qJ(1));
t26 = g(1) * t23 - g(2) * t24;
t25 = -t15 * t10 - t16 * t22;
t18 = pkin(7) + qJ(4);
t17 = t24 * pkin(1);
t14 = cos(t18);
t13 = sin(t18);
t4 = t5 * t21;
t3 = t5 * sin(pkin(7));
t2 = t5 * t14;
t1 = t5 * t13;
t7 = [0, 0, 0, 0, 0, 0, t26, g(1) * t24 + g(2) * t23, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t26 * pkin(1), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (t28 - t30) - g(2) * (t17 + t29), 0, 0, 0, 0, 0, 0, t2, -t1, -t6, -g(1) * (t25 - t30) - g(2) * (t17 + t27); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * t28 - g(2) * t29, 0, 0, 0, 0, 0, 0, t2, -t1, -t6, -g(1) * t25 - g(2) * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t14 + t6 * t13, g(3) * t13 + t6 * t14, 0, 0;];
taug_reg = t7;
