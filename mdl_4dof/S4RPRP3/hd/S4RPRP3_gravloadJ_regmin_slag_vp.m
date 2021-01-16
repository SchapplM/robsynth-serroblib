% Calculate minimal parameter regressor of gravitation load for
% S4RPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% taug_reg [4x15]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:20
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RPRP3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:20:36
% EndTime: 2021-01-15 10:20:37
% DurationCPUTime: 0.09s
% Computational Cost: add. (56->20), mult. (66->25), div. (0->0), fcn. (59->6), ass. (0->17)
t8 = qJ(1) + pkin(6);
t6 = sin(t8);
t7 = cos(t8);
t16 = g(1) * t7 + g(2) * t6;
t15 = g(1) * t6 - g(2) * t7;
t11 = sin(qJ(1));
t13 = cos(qJ(1));
t14 = g(1) * t11 - g(2) * t13;
t10 = sin(qJ(3));
t12 = cos(qJ(3));
t1 = -g(3) * t12 + t16 * t10;
t9 = -qJ(4) - pkin(5);
t5 = t12 * pkin(3) + pkin(2);
t4 = t15 * t12;
t3 = t15 * t10;
t2 = g(3) * t10 + t16 * t12;
t17 = [0, t14, g(1) * t13 + g(2) * t11, t14 * pkin(1), 0, 0, 0, 0, 0, t4, -t3, t4, -t3, -t16, -g(1) * (-t11 * pkin(1) - t6 * t5 - t7 * t9) - g(2) * (t13 * pkin(1) + t7 * t5 - t6 * t9); 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, t2, 0, t1 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15;];
taug_reg = t17;
