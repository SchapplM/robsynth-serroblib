% Calculate minimal parameter regressor of gravitation load for
% S5RRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% taug_reg [5x28]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR12_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR12_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:30:37
% EndTime: 2019-12-31 20:30:37
% DurationCPUTime: 0.14s
% Computational Cost: add. (99->37), mult. (274->60), div. (0->0), fcn. (315->8), ass. (0->32)
t16 = sin(qJ(5));
t17 = sin(qJ(4));
t18 = sin(qJ(2));
t21 = cos(qJ(2));
t33 = cos(qJ(4));
t11 = t18 * t17 + t21 * t33;
t19 = sin(qJ(1));
t37 = -t21 * t17 + t18 * t33;
t5 = t37 * t19;
t22 = cos(qJ(1));
t7 = t37 * t22;
t24 = -g(1) * t7 - g(2) * t5 + g(3) * t11;
t39 = t24 * t16;
t20 = cos(qJ(5));
t38 = t24 * t20;
t13 = g(1) * t22 + g(2) * t19;
t34 = g(3) * t37;
t30 = g(1) * t19 - g(2) * t22;
t6 = t11 * t19;
t29 = t6 * t16 - t22 * t20;
t28 = t22 * t16 + t6 * t20;
t27 = t21 * pkin(2) + t18 * qJ(3);
t25 = pkin(1) + t27;
t8 = t11 * t22;
t23 = g(1) * t8 + g(2) * t6 + t34;
t10 = t30 * t21;
t9 = t30 * t18;
t4 = g(3) * t18 + t13 * t21;
t3 = -g(3) * t21 + t13 * t18;
t2 = -t19 * t16 + t8 * t20;
t1 = -t8 * t16 - t19 * t20;
t12 = [0, t30, t13, 0, 0, 0, 0, 0, t10, -t9, t10, -t13, t9, (-g(1) * pkin(6) - g(2) * t25) * t22 + (-g(2) * pkin(6) + g(1) * t25) * t19, 0, 0, 0, 0, 0, g(1) * t6 - g(2) * t8, g(1) * t5 - g(2) * t7, 0, 0, 0, 0, 0, g(1) * t28 - g(2) * t2, -g(1) * t29 - g(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, t3, 0, -t4, -g(3) * t27 + t13 * (pkin(2) * t18 - qJ(3) * t21), 0, 0, 0, 0, 0, -t24, -t23, 0, 0, 0, 0, 0, -t38, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t23, 0, 0, 0, 0, 0, t38, -t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 + g(2) * t29 + t16 * t34, g(1) * t2 + g(2) * t28 + t20 * t34;];
taug_reg = t12;
