% Calculate minimal parameter regressor of gravitation load for
% S5PRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% taug_reg [5x22]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:24
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRRP3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:23:27
% EndTime: 2021-01-15 16:23:27
% DurationCPUTime: 0.09s
% Computational Cost: add. (138->23), mult. (111->29), div. (0->0), fcn. (103->6), ass. (0->18)
t17 = qJ(3) + qJ(4);
t13 = cos(t17);
t19 = cos(qJ(3));
t20 = t19 * pkin(3) + pkin(4) * t13;
t16 = pkin(8) + qJ(2);
t10 = sin(t16);
t11 = cos(t16);
t6 = g(1) * t11 + g(2) * t10;
t5 = g(1) * t10 - g(2) * t11;
t12 = sin(t17);
t1 = -g(3) * t13 + t6 * t12;
t18 = sin(qJ(3));
t15 = -qJ(5) - pkin(7) - pkin(6);
t7 = pkin(2) + t20;
t4 = t5 * t13;
t3 = t5 * t12;
t2 = g(3) * t12 + t6 * t13;
t8 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, t5, t6, 0, 0, 0, 0, 0, t5 * t19, -t5 * t18, 0, 0, 0, 0, 0, t4, -t3, t4, -t3, -t6, -g(1) * (-t10 * t7 - t11 * t15) - g(2) * (-t10 * t15 + t11 * t7); 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t19 + t6 * t18, g(3) * t18 + t6 * t19, 0, 0, 0, 0, 0, t1, t2, t1, t2, 0, -g(3) * t20 - t6 * (-t18 * pkin(3) - pkin(4) * t12); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5;];
taug_reg = t8;
