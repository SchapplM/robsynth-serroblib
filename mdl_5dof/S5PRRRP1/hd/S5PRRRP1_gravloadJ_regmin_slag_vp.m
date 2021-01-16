% Calculate minimal parameter regressor of gravitation load for
% S5PRRRP1
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
% taug_reg [5x18]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:15
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRRP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:14:40
% EndTime: 2021-01-15 16:14:41
% DurationCPUTime: 0.08s
% Computational Cost: add. (147->25), mult. (96->26), div. (0->0), fcn. (89->6), ass. (0->19)
t15 = pkin(8) + qJ(2);
t14 = qJ(3) + t15;
t10 = cos(t14);
t18 = cos(qJ(4));
t11 = t18 * pkin(4) + pkin(3);
t16 = qJ(5) + pkin(7);
t9 = sin(t14);
t20 = t16 * t10 - t9 * t11;
t19 = t10 * t11 + t9 * t16;
t5 = g(1) * t9 - g(2) * t10;
t6 = g(1) * t10 + g(2) * t9;
t17 = sin(qJ(4));
t1 = -g(3) * t18 + t6 * t17;
t13 = cos(t15);
t12 = sin(t15);
t4 = t5 * t18;
t3 = t5 * t17;
t2 = g(3) * t17 + t6 * t18;
t7 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, g(1) * t12 - g(2) * t13, g(1) * t13 + g(2) * t12, 0, t5, t6, 0, 0, 0, 0, 0, t4, -t3, t4, -t3, -t6, -g(1) * (-pkin(2) * t12 + t20) - g(2) * (pkin(2) * t13 + t19); 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t4, -t3, t4, -t3, -t6, -g(1) * t20 - g(2) * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5;];
taug_reg = t7;
