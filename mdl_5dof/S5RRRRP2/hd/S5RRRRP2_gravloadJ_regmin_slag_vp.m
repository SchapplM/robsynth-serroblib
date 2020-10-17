% Calculate minimal parameter regressor of gravitation load for
% S5RRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% taug_reg [5x22]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRP2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:11:45
% EndTime: 2020-01-03 12:11:46
% DurationCPUTime: 0.09s
% Computational Cost: add. (144->29), mult. (121->37), div. (0->0), fcn. (111->8), ass. (0->24)
t22 = qJ(3) + qJ(4);
t18 = cos(t22);
t26 = cos(qJ(3));
t29 = t26 * pkin(3) + pkin(4) * t18;
t11 = pkin(2) + t29;
t23 = qJ(1) + qJ(2);
t17 = sin(t23);
t19 = cos(t23);
t21 = -qJ(5) - pkin(8) - pkin(7);
t30 = t17 * t11 + t19 * t21;
t28 = t19 * t11 - t17 * t21;
t10 = g(2) * t19 + g(3) * t17;
t9 = g(2) * t17 - g(3) * t19;
t16 = sin(t22);
t2 = -g(1) * t18 + t9 * t16;
t27 = cos(qJ(1));
t25 = sin(qJ(1));
t24 = sin(qJ(3));
t6 = t10 * t26;
t5 = t10 * t24;
t4 = t10 * t18;
t3 = t10 * t16;
t1 = g(1) * t16 + t9 * t18;
t7 = [0, -g(2) * t27 - g(3) * t25, g(2) * t25 - g(3) * t27, 0, -t10, t9, 0, 0, 0, 0, 0, -t6, t5, 0, 0, 0, 0, 0, -t4, t3, -t9, -g(2) * (t27 * pkin(1) + t28) - g(3) * (t25 * pkin(1) + t30); 0, 0, 0, 0, -t10, t9, 0, 0, 0, 0, 0, -t6, t5, 0, 0, 0, 0, 0, -t4, t3, -t9, -g(2) * t28 - g(3) * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t26 + t9 * t24, g(1) * t24 + t9 * t26, 0, 0, 0, 0, 0, t2, t1, 0, -g(1) * t29 - t9 * (-t24 * pkin(3) - pkin(4) * t16); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10;];
taug_reg = t7;
