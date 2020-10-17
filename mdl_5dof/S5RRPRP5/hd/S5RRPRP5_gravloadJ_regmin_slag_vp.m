% Calculate minimal parameter regressor of gravitation load for
% S5RRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% taug_reg [5x23]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRP5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:54:56
% EndTime: 2019-12-31 19:54:57
% DurationCPUTime: 0.13s
% Computational Cost: add. (181->40), mult. (159->49), div. (0->0), fcn. (146->8), ass. (0->28)
t25 = cos(qJ(2));
t18 = t25 * pkin(2);
t21 = qJ(2) + pkin(8);
t17 = qJ(4) + t21;
t14 = sin(t17);
t15 = cos(t17);
t32 = t15 * pkin(4) + t14 * qJ(5);
t34 = pkin(3) * cos(t21) + t18 + t32;
t33 = pkin(4) * t14;
t22 = -qJ(3) - pkin(6);
t30 = qJ(5) * t15;
t23 = sin(qJ(2));
t29 = -pkin(3) * sin(t21) - t23 * pkin(2) - t33;
t24 = sin(qJ(1));
t26 = cos(qJ(1));
t10 = g(1) * t26 + g(2) * t24;
t9 = g(1) * t24 - g(2) * t26;
t28 = pkin(1) + t34;
t27 = -g(3) * t25 + t10 * t23;
t20 = -pkin(7) + t22;
t16 = t18 + pkin(1);
t8 = t26 * t30;
t7 = t24 * t30;
t4 = t9 * t15;
t3 = t9 * t14;
t2 = g(3) * t14 + t10 * t15;
t1 = -g(3) * t15 + t10 * t14;
t5 = [0, t9, t10, 0, 0, 0, 0, 0, t9 * t25, -t9 * t23, -t10, -g(1) * (-t24 * t16 - t26 * t22) - g(2) * (t26 * t16 - t24 * t22), 0, 0, 0, 0, 0, t4, -t3, t4, -t10, t3, (g(1) * t20 - g(2) * t28) * t26 + (g(1) * t28 + g(2) * t20) * t24; 0, 0, 0, 0, 0, 0, 0, 0, t27, g(3) * t23 + t10 * t25, 0, t27 * pkin(2), 0, 0, 0, 0, 0, t1, t2, t1, 0, -t2, -g(1) * (t29 * t26 + t8) - g(2) * (t29 * t24 + t7) - g(3) * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, 0, -t2, -g(1) * (-t26 * t33 + t8) - g(2) * (-t24 * t33 + t7) - g(3) * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t5;
