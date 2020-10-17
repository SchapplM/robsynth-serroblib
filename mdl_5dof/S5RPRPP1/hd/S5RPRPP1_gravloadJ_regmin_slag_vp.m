% Calculate minimal parameter regressor of gravitation load for
% S5RPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
% 
% Output:
% taug_reg [5x17]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:09
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:09:14
% EndTime: 2019-12-31 18:09:14
% DurationCPUTime: 0.11s
% Computational Cost: add. (125->36), mult. (111->44), div. (0->0), fcn. (99->8), ass. (0->22)
t12 = qJ(3) + pkin(8);
t6 = sin(t12);
t8 = cos(t12);
t22 = t8 * pkin(4) + t6 * qJ(5);
t13 = qJ(1) + pkin(7);
t7 = sin(t13);
t9 = cos(t13);
t25 = g(1) * t9 + g(2) * t7;
t18 = cos(qJ(1));
t17 = cos(qJ(3));
t10 = t17 * pkin(3);
t5 = t10 + pkin(2);
t27 = t18 * pkin(1) + t9 * t5;
t24 = g(1) * t7 - g(2) * t9;
t14 = -qJ(4) - pkin(6);
t16 = sin(qJ(1));
t23 = -t16 * pkin(1) - t9 * t14;
t21 = g(1) * t16 - g(2) * t18;
t15 = sin(qJ(3));
t19 = -g(3) * t17 + t25 * t15;
t1 = -g(3) * t8 + t25 * t6;
t2 = [0, t21, g(1) * t18 + g(2) * t16, t21 * pkin(1), 0, 0, 0, 0, 0, t24 * t17, -t24 * t15, -t25, -g(1) * (-t7 * t5 + t23) - g(2) * (-t7 * t14 + t27), t24 * t8, -t25, t24 * t6, -g(1) * t23 - g(2) * (t22 * t9 + t27) + (-g(1) * (-t22 - t5) + g(2) * t14) * t7; 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, g(3) * t15 + t25 * t17, 0, t19 * pkin(3), t1, 0, -g(3) * t6 - t25 * t8, -g(3) * (t10 + t22) + t25 * (pkin(3) * t15 + pkin(4) * t6 - qJ(5) * t8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, 0, 0, 0, -t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t2;
