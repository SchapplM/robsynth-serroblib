% Calculate minimal parameter regressor of gravitation load for
% S5RPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% taug_reg [5x24]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:17
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR14_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR14_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:17:03
% EndTime: 2021-01-15 12:17:04
% DurationCPUTime: 0.12s
% Computational Cost: add. (78->36), mult. (126->50), div. (0->0), fcn. (127->8), ass. (0->25)
t15 = sin(qJ(1));
t18 = cos(qJ(1));
t5 = g(1) * t15 - g(2) * t18;
t12 = qJ(3) + pkin(8);
t8 = sin(t12);
t9 = cos(t12);
t27 = -g(3) * t8 + t5 * t9;
t25 = g(3) * t9;
t13 = sin(qJ(5));
t24 = t15 * t13;
t16 = cos(qJ(5));
t23 = t15 * t16;
t22 = t18 * t13;
t21 = t18 * t16;
t6 = g(1) * t18 + g(2) * t15;
t14 = sin(qJ(3));
t17 = cos(qJ(3));
t19 = g(3) * t14 - t5 * t17;
t11 = pkin(1) + pkin(6) + qJ(4);
t7 = t14 * pkin(3) + qJ(2);
t4 = t8 * t21 - t24;
t3 = t8 * t22 + t23;
t2 = t8 * t23 + t22;
t1 = -t8 * t24 + t21;
t10 = [0, t5, t6, -t5, -t6, -g(1) * (-t15 * pkin(1) + t18 * qJ(2)) - g(2) * (t18 * pkin(1) + t15 * qJ(2)), 0, 0, 0, 0, 0, -t6 * t14, -t6 * t17, -t6 * t8, -t6 * t9, t5, -g(1) * (-t11 * t15 + t7 * t18) - g(2) * (t11 * t18 + t7 * t15), 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t2, g(1) * t3 - g(2) * t1; 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, g(3) * t17 + t5 * t14, -t27, t5 * t8 + t25, 0, t19 * pkin(3), 0, 0, 0, 0, 0, -t27 * t16, t27 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 - g(2) * t3 + t13 * t25, g(1) * t2 - g(2) * t4 + t16 * t25;];
taug_reg = t10;
