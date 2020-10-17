% Calculate minimal parameter regressor of gravitation load for
% S5RPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% taug_reg [5x22]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:57
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRP12_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP12_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:57:20
% EndTime: 2019-12-31 18:57:21
% DurationCPUTime: 0.15s
% Computational Cost: add. (68->36), mult. (155->52), div. (0->0), fcn. (153->6), ass. (0->28)
t18 = sin(qJ(1));
t21 = cos(qJ(1));
t38 = -g(1) * t18 + g(2) * t21;
t16 = sin(qJ(4));
t17 = sin(qJ(3));
t19 = cos(qJ(4));
t27 = t21 * t19;
t30 = t18 * t16;
t3 = -t17 * t30 + t27;
t20 = cos(qJ(3));
t31 = g(3) * t20;
t28 = t21 * t16;
t29 = t18 * t19;
t5 = t17 * t28 + t29;
t37 = -g(1) * t3 - g(2) * t5 + t16 * t31;
t2 = -g(3) * t17 - t20 * t38;
t25 = pkin(4) * t16 + pkin(6);
t24 = g(2) * (t21 * pkin(1) + t18 * qJ(2));
t9 = g(1) * t21 + g(2) * t18;
t10 = t19 * pkin(4) + pkin(3);
t15 = -qJ(5) - pkin(7);
t22 = t17 * t10 + t20 * t15;
t12 = t21 * qJ(2);
t7 = t9 * t20;
t6 = t17 * t27 - t30;
t4 = t17 * t29 + t28;
t1 = -t17 * t38 + t31;
t8 = [0, -t38, t9, t38, -t9, -g(1) * (-t18 * pkin(1) + t12) - t24, 0, 0, 0, 0, 0, -t9 * t17, -t7, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t4, g(1) * t5 - g(2) * t3, t7, -g(1) * t12 - t24 + (-g(1) * t22 - g(2) * t25) * t21 + (-g(1) * (-pkin(1) - t25) - g(2) * t22) * t18; 0, 0, 0, 0, 0, t38, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, 0, 0, 0, 0, 0, -t2 * t19, t2 * t16, -t1, g(3) * t22 + t38 * (t10 * t20 - t15 * t17); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, g(1) * t4 - g(2) * t6 + t19 * t31, 0, t37 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2;];
taug_reg = t8;
