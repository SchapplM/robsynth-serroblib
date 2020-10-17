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
% taug_reg [5x20]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 16:44:13
% EndTime: 2019-12-05 16:44:13
% DurationCPUTime: 0.08s
% Computational Cost: add. (100->22), mult. (83->29), div. (0->0), fcn. (75->6), ass. (0->16)
t15 = qJ(3) + qJ(4);
t11 = cos(t15);
t17 = cos(qJ(3));
t18 = t17 * pkin(3) + pkin(4) * t11;
t14 = pkin(8) + qJ(2);
t8 = sin(t14);
t9 = cos(t14);
t4 = g(1) * t9 + g(2) * t8;
t3 = g(1) * t8 - g(2) * t9;
t10 = sin(t15);
t1 = -g(3) * t11 + t4 * t10;
t16 = sin(qJ(3));
t13 = -qJ(5) - pkin(7) - pkin(6);
t5 = pkin(2) + t18;
t2 = g(3) * t10 + t4 * t11;
t6 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, t3, t4, 0, 0, 0, 0, 0, t3 * t17, -t3 * t16, 0, 0, 0, 0, 0, t3 * t11, -t3 * t10, -t4, -g(1) * (-t9 * t13 - t8 * t5) - g(2) * (-t8 * t13 + t9 * t5); 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t17 + t4 * t16, g(3) * t16 + t4 * t17, 0, 0, 0, 0, 0, t1, t2, 0, -g(3) * t18 - t4 * (-t16 * pkin(3) - pkin(4) * t10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3;];
taug_reg = t6;
