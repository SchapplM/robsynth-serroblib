% Calculate minimal parameter regressor of gravitation load for
% S5RRRRP1
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
% taug_reg [5x28]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:53
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 23:52:43
% EndTime: 2021-01-15 23:52:43
% DurationCPUTime: 0.11s
% Computational Cost: add. (212->27), mult. (175->39), div. (0->0), fcn. (164->8), ass. (0->24)
t21 = qJ(2) + qJ(3);
t19 = qJ(4) + t21;
t15 = cos(t19);
t17 = cos(t21);
t28 = pkin(3) * t17 + pkin(4) * t15;
t24 = cos(qJ(2));
t27 = t24 * pkin(2) + t28;
t14 = sin(t19);
t16 = sin(t21);
t26 = -pkin(3) * t16 - pkin(4) * t14;
t23 = sin(qJ(1));
t25 = cos(qJ(1));
t11 = g(1) * t25 + g(2) * t23;
t10 = g(1) * t23 - g(2) * t25;
t1 = -g(3) * t15 + t11 * t14;
t22 = sin(qJ(2));
t18 = -qJ(5) - pkin(8) - pkin(7) - pkin(6);
t7 = pkin(1) + t27;
t6 = t10 * t15;
t5 = t10 * t14;
t4 = g(3) * t16 + t11 * t17;
t3 = -g(3) * t17 + t11 * t16;
t2 = g(3) * t14 + t11 * t15;
t8 = [0, t10, t11, 0, 0, 0, 0, 0, t10 * t24, -t10 * t22, 0, 0, 0, 0, 0, t10 * t17, -t10 * t16, 0, 0, 0, 0, 0, t6, -t5, t6, -t5, -t11, -g(1) * (-t25 * t18 - t23 * t7) - g(2) * (-t23 * t18 + t25 * t7); 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t24 + t11 * t22, g(3) * t22 + t11 * t24, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2, t1, t2, 0, -g(3) * t27 - t11 * (-t22 * pkin(2) + t26); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2, t1, t2, 0, -g(3) * t28 - t11 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10;];
taug_reg = t8;
