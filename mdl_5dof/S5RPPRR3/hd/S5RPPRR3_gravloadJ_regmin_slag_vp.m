% Calculate minimal parameter regressor of gravitation load for
% S5RPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% taug_reg [5x22]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRR3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:28:14
% EndTime: 2020-01-03 11:28:15
% DurationCPUTime: 0.10s
% Computational Cost: add. (104->22), mult. (74->31), div. (0->0), fcn. (70->10), ass. (0->17)
t10 = pkin(9) + qJ(4);
t11 = qJ(1) + pkin(8);
t6 = sin(t11);
t8 = cos(t11);
t18 = g(2) * t8 + g(3) * t6;
t17 = g(2) * t6 - g(3) * t8;
t14 = sin(qJ(1));
t15 = cos(qJ(1));
t16 = -g(2) * t15 - g(3) * t14;
t9 = qJ(5) + t10;
t7 = cos(t10);
t5 = sin(t10);
t4 = cos(t9);
t3 = sin(t9);
t2 = -g(1) * t4 + t17 * t3;
t1 = g(1) * t3 + t17 * t4;
t12 = [0, t16, g(2) * t14 - g(3) * t15, t16 * pkin(1), -t18 * cos(pkin(9)), t18 * sin(pkin(9)), -t17, -g(2) * (t15 * pkin(1) + t8 * pkin(2) + t6 * qJ(3)) - g(3) * (t14 * pkin(1) + t6 * pkin(2) - t8 * qJ(3)), 0, 0, 0, 0, 0, -t18 * t7, t18 * t5, 0, 0, 0, 0, 0, -t18 * t4, t18 * t3; 0, 0, 0, -g(1), 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 + t17 * t5, g(1) * t5 + t17 * t7, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1;];
taug_reg = t12;
