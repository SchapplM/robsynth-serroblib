% Calculate minimal parameter regressor of gravitation load for
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% taug_reg [5x17]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:07:36
% EndTime: 2020-01-03 12:07:37
% DurationCPUTime: 0.08s
% Computational Cost: add. (128->23), mult. (76->33), div. (0->0), fcn. (70->10), ass. (0->24)
t19 = qJ(1) + qJ(2);
t18 = qJ(3) + t19;
t14 = sin(t18);
t16 = sin(t19);
t27 = pkin(2) * t16 + pkin(3) * t14;
t15 = cos(t18);
t17 = cos(t19);
t26 = pkin(2) * t17 + pkin(3) * t15;
t13 = pkin(9) + t18;
t7 = sin(t13);
t8 = cos(t13);
t25 = g(2) * t8 + g(3) * t7;
t24 = g(2) * t7 - g(3) * t8;
t4 = -g(2) * t15 - g(3) * t14;
t23 = cos(qJ(1));
t22 = cos(qJ(5));
t21 = sin(qJ(1));
t20 = sin(qJ(5));
t6 = -g(2) * t17 - g(3) * t16;
t5 = g(2) * t16 - g(3) * t17;
t3 = g(2) * t14 - g(3) * t15;
t2 = t25 * t22;
t1 = t25 * t20;
t9 = [0, -g(2) * t23 - g(3) * t21, g(2) * t21 - g(3) * t23, 0, t6, t5, 0, t4, t3, -g(2) * (t23 * pkin(1) + t26) - g(3) * (t21 * pkin(1) + t27), 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, t6, t5, 0, t4, t3, -g(2) * t26 - g(3) * t27, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, t4, t3, t4 * pkin(3), 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t22 + t24 * t20, g(1) * t20 + t24 * t22;];
taug_reg = t9;
