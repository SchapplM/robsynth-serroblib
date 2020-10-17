% Calculate minimal parameter regressor of gravitation load for
% S5RRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% taug_reg [5x19]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPPR4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:27:49
% EndTime: 2019-12-31 19:27:49
% DurationCPUTime: 0.11s
% Computational Cost: add. (144->34), mult. (130->38), div. (0->0), fcn. (142->8), ass. (0->25)
t21 = sin(qJ(1));
t31 = t21 * pkin(1);
t18 = qJ(1) + qJ(2);
t15 = sin(t18);
t16 = cos(t18);
t30 = t16 * pkin(2) + t15 * qJ(3);
t29 = cos(pkin(8));
t23 = cos(qJ(1));
t28 = t23 * pkin(1) + t30;
t11 = t16 * qJ(3);
t27 = -t15 * pkin(2) + t11;
t19 = sin(pkin(8));
t5 = -t15 * t29 + t16 * t19;
t6 = t15 * t19 + t16 * t29;
t26 = g(1) * t6 - g(2) * t5;
t25 = g(1) * t5 + g(2) * t6;
t24 = t11 + (-pkin(2) - pkin(3)) * t15;
t22 = cos(qJ(5));
t20 = sin(qJ(5));
t12 = t16 * pkin(3);
t8 = g(1) * t16 + g(2) * t15;
t7 = g(1) * t15 - g(2) * t16;
t2 = t25 * t22;
t1 = t25 * t20;
t3 = [0, g(1) * t21 - g(2) * t23, g(1) * t23 + g(2) * t21, 0, t7, t8, t7, -t8, -g(1) * (t27 - t31) - g(2) * t28, -t25, -t26, -g(1) * (t24 - t31) - g(2) * (t12 + t28), 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, t7, t8, t7, -t8, -g(1) * t27 - g(2) * t30, -t25, -t26, -g(1) * t24 - g(2) * (t12 + t30), 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t22 + t26 * t20, -g(3) * t20 + t26 * t22;];
taug_reg = t3;
