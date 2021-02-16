% Calculate minimal parameter regressor of gravitation load for
% S5RPRPR1
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
% Datum: 2021-01-15 11:34
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:33:52
% EndTime: 2021-01-15 11:33:52
% DurationCPUTime: 0.09s
% Computational Cost: add. (86->29), mult. (102->36), div. (0->0), fcn. (95->8), ass. (0->18)
t13 = qJ(3) + pkin(8);
t15 = sin(qJ(1));
t17 = cos(qJ(1));
t4 = g(1) * t17 + g(2) * t15;
t3 = g(1) * t15 - g(2) * t17;
t14 = sin(qJ(3));
t16 = cos(qJ(3));
t18 = g(3) * t14 - t3 * t16;
t12 = pkin(1) + pkin(6) + qJ(4);
t10 = qJ(5) + t13;
t9 = cos(t13);
t8 = sin(t13);
t7 = t14 * pkin(3) + qJ(2);
t6 = cos(t10);
t5 = sin(t10);
t2 = g(3) * t5 - t3 * t6;
t1 = g(3) * t6 + t3 * t5;
t11 = [0, t3, t4, -t3, -t4, -g(1) * (-t15 * pkin(1) + t17 * qJ(2)) - g(2) * (t17 * pkin(1) + t15 * qJ(2)), 0, 0, 0, 0, 0, -t4 * t14, -t4 * t16, -t4 * t8, -t4 * t9, t3, -g(1) * (-t12 * t15 + t7 * t17) - g(2) * (t12 * t17 + t7 * t15), 0, 0, 0, 0, 0, -t4 * t5, -t4 * t6; 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, g(3) * t16 + t3 * t14, g(3) * t8 - t3 * t9, g(3) * t9 + t3 * t8, 0, t18 * pkin(3), 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1;];
taug_reg = t11;
