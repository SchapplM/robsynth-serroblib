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
% taug_reg [5x24]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 00:01
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
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
% StartTime: 2021-01-16 00:01:01
% EndTime: 2021-01-16 00:01:02
% DurationCPUTime: 0.11s
% Computational Cost: add. (192->31), mult. (157->37), div. (0->0), fcn. (147->8), ass. (0->24)
t21 = qJ(3) + qJ(4);
t17 = cos(t21);
t25 = cos(qJ(3));
t29 = t25 * pkin(3) + pkin(4) * t17;
t10 = -pkin(2) - t29;
t22 = qJ(1) + qJ(2);
t16 = sin(t22);
t18 = cos(t22);
t20 = -qJ(5) - pkin(8) - pkin(7);
t28 = -t18 * t10 - t16 * t20;
t27 = -t10 * t16 + t18 * t20;
t9 = g(2) * t18 + g(3) * t16;
t8 = g(2) * t16 - g(3) * t18;
t15 = sin(t21);
t2 = -g(1) * t17 + t8 * t15;
t26 = cos(qJ(1));
t24 = sin(qJ(1));
t23 = sin(qJ(3));
t6 = t9 * t25;
t5 = t9 * t23;
t4 = t9 * t17;
t3 = t9 * t15;
t1 = g(1) * t15 + t8 * t17;
t7 = [0, -g(2) * t26 - g(3) * t24, g(2) * t24 - g(3) * t26, 0, -t9, t8, 0, 0, 0, 0, 0, -t6, t5, 0, 0, 0, 0, 0, -t4, t3, -t4, t3, -t8, -g(2) * (t26 * pkin(1) + t28) - g(3) * (t24 * pkin(1) + t27); 0, 0, 0, 0, -t9, t8, 0, 0, 0, 0, 0, -t6, t5, 0, 0, 0, 0, 0, -t4, t3, -t4, t3, -t8, -g(2) * t28 - g(3) * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t25 + t8 * t23, g(1) * t23 + t8 * t25, 0, 0, 0, 0, 0, t2, t1, t2, t1, 0, -g(1) * t29 + t8 * (t23 * pkin(3) + pkin(4) * t15); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, t2, t1, 0, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9;];
taug_reg = t7;
