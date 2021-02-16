% Calculate minimal parameter regressor of gravitation load for
% S4RRPR6
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
% taug_reg [4x21]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:46
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRPR6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:46:22
% EndTime: 2021-01-15 10:46:22
% DurationCPUTime: 0.08s
% Computational Cost: add. (77->20), mult. (88->30), div. (0->0), fcn. (83->8), ass. (0->18)
t11 = qJ(2) + pkin(7);
t14 = sin(qJ(1));
t16 = cos(qJ(1));
t4 = g(1) * t16 + g(2) * t14;
t3 = g(1) * t14 - g(2) * t16;
t13 = sin(qJ(2));
t15 = cos(qJ(2));
t17 = -g(3) * t15 + t4 * t13;
t12 = -qJ(3) - pkin(5);
t10 = qJ(4) + t11;
t9 = cos(t11);
t8 = sin(t11);
t7 = t15 * pkin(2) + pkin(1);
t6 = cos(t10);
t5 = sin(t10);
t2 = g(3) * t5 + t4 * t6;
t1 = -g(3) * t6 + t4 * t5;
t18 = [0, t3, t4, 0, 0, 0, 0, 0, t3 * t15, -t3 * t13, t3 * t9, -t3 * t8, -t4, -g(1) * (-t12 * t16 - t14 * t7) - g(2) * (-t14 * t12 + t16 * t7), 0, 0, 0, 0, 0, t3 * t6, -t3 * t5; 0, 0, 0, 0, 0, 0, 0, 0, t17, g(3) * t13 + t4 * t15, -g(3) * t9 + t4 * t8, g(3) * t8 + t4 * t9, 0, t17 * pkin(2), 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t18;
