% Calculate minimal parameter regressor of gravitation load for
% S4RRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% 
% Output:
% taug_reg [4x18]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:36
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRPP3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:36:05
% EndTime: 2021-01-15 10:36:06
% DurationCPUTime: 0.09s
% Computational Cost: add. (87->33), mult. (129->43), div. (0->0), fcn. (121->8), ass. (0->26)
t18 = qJ(2) + pkin(6);
t15 = cos(t18);
t28 = g(3) * t15;
t21 = -qJ(3) - pkin(5);
t23 = sin(qJ(1));
t27 = t23 * t21;
t25 = cos(qJ(1));
t11 = g(1) * t25 + g(2) * t23;
t10 = g(1) * t23 - g(2) * t25;
t22 = sin(qJ(2));
t24 = cos(qJ(2));
t26 = -g(3) * t24 + t11 * t22;
t20 = cos(pkin(6));
t19 = sin(pkin(6));
t16 = t24 * pkin(2);
t14 = sin(t18);
t13 = t16 + pkin(1);
t12 = t21 * t25;
t9 = -t19 * pkin(3) + qJ(4) * t20;
t8 = pkin(3) * t20 + qJ(4) * t19 + pkin(2);
t6 = t10 * t15;
t5 = t10 * t14;
t4 = g(3) * t14 + t11 * t15;
t3 = t11 * t14 - t28;
t1 = t9 * t22 + t8 * t24 + pkin(1);
t2 = [0, t10, t11, 0, 0, 0, 0, 0, t10 * t24, -t10 * t22, t6, -t5, -t11, -g(1) * (-t23 * t13 - t12) - g(2) * (t25 * t13 - t27), t6, -t11, t5, -g(1) * (-t1 * t23 - t12) - g(2) * (t1 * t25 - t27); 0, 0, 0, 0, 0, 0, 0, 0, t26, g(3) * t22 + t11 * t24, t3, t4, 0, t26 * pkin(2), t3, 0, -t4, -g(3) * (t15 * pkin(3) + t14 * qJ(4) + t16) - t11 * (-t8 * t22 + t9 * t24); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28 - t11 * (t19 * t24 + t20 * t22);];
taug_reg = t2;
