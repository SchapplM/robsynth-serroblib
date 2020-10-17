% Calculate minimal parameter regressor of gravitation load for
% S5RPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% taug_reg [5x18]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:34:01
% EndTime: 2020-01-03 11:34:02
% DurationCPUTime: 0.09s
% Computational Cost: add. (134->27), mult. (80->29), div. (0->0), fcn. (74->10), ass. (0->19)
t17 = qJ(1) + pkin(8);
t15 = qJ(3) + t17;
t11 = sin(t15);
t12 = cos(t15);
t24 = t12 * pkin(3) + t11 * qJ(4);
t23 = t11 * pkin(3) - t12 * qJ(4);
t6 = g(2) * t12 + g(3) * t11;
t5 = g(2) * t11 - g(3) * t12;
t20 = sin(qJ(1));
t21 = cos(qJ(1));
t22 = -g(2) * t21 - g(3) * t20;
t16 = pkin(9) + qJ(5);
t14 = cos(t16);
t13 = sin(t16);
t4 = t6 * cos(pkin(9));
t3 = t6 * sin(pkin(9));
t2 = t6 * t14;
t1 = t6 * t13;
t7 = [0, t22, g(2) * t20 - g(3) * t21, t22 * pkin(1), 0, -t6, t5, -t4, t3, -t5, -g(2) * (pkin(2) * cos(t17) + t21 * pkin(1) + t24) - g(3) * (pkin(2) * sin(t17) + t20 * pkin(1) + t23), 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t6, t5, -t4, t3, -t5, -g(2) * t24 - g(3) * t23, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t14 + t5 * t13, g(1) * t13 + t5 * t14;];
taug_reg = t7;
