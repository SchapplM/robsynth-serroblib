% Calculate minimal parameter regressor of gravitation load for
% S4RPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% 
% Output:
% taug_reg [4x17]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:27
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RPRP6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_gravloadJ_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:27:40
% EndTime: 2021-01-15 10:27:41
% DurationCPUTime: 0.07s
% Computational Cost: add. (34->21), mult. (74->24), div. (0->0), fcn. (67->4), ass. (0->13)
t11 = sin(qJ(1));
t13 = cos(qJ(1));
t6 = g(1) * t13 + g(2) * t11;
t5 = g(1) * t11 - g(2) * t13;
t10 = sin(qJ(3));
t12 = cos(qJ(3));
t2 = g(3) * t10 - t5 * t12;
t9 = pkin(1) + pkin(5) + qJ(4);
t7 = t10 * pkin(3) + qJ(2);
t4 = t6 * t12;
t3 = t6 * t10;
t1 = g(3) * t12 + t5 * t10;
t8 = [0, t5, t6, -t5, -t6, -g(1) * (-t11 * pkin(1) + t13 * qJ(2)) - g(2) * (t13 * pkin(1) + t11 * qJ(2)), 0, 0, 0, 0, 0, -t3, -t4, -t3, -t4, t5, -g(1) * (-t9 * t11 + t7 * t13) - g(2) * (t7 * t11 + t9 * t13); 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, t2, t1, 0, t2 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6;];
taug_reg = t8;
