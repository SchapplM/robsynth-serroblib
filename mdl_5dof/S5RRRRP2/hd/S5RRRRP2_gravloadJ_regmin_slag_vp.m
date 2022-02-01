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
% Datum: 2022-01-20 11:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 11:49:27
% EndTime: 2022-01-20 11:49:28
% DurationCPUTime: 0.10s
% Computational Cost: add. (192->30), mult. (157->37), div. (0->0), fcn. (147->8), ass. (0->24)
t20 = qJ(3) + qJ(4);
t16 = cos(t20);
t24 = cos(qJ(3));
t28 = t24 * pkin(3) + pkin(4) * t16;
t10 = pkin(2) + t28;
t21 = qJ(1) + qJ(2);
t15 = sin(t21);
t17 = cos(t21);
t19 = qJ(5) + pkin(7) + pkin(8);
t27 = t17 * t10 + t15 * t19;
t26 = -t15 * t10 + t19 * t17;
t9 = g(1) * t17 + g(2) * t15;
t8 = g(1) * t15 - g(2) * t17;
t14 = sin(t20);
t1 = -g(3) * t16 + t9 * t14;
t25 = cos(qJ(1));
t23 = sin(qJ(1));
t22 = sin(qJ(3));
t6 = t8 * t24;
t5 = t8 * t22;
t4 = t8 * t16;
t3 = t8 * t14;
t2 = g(3) * t14 + t9 * t16;
t7 = [0, g(1) * t23 - g(2) * t25, g(1) * t25 + g(2) * t23, 0, t8, t9, 0, 0, 0, 0, 0, t6, -t5, 0, 0, 0, 0, 0, t4, -t3, t4, -t3, -t9, -g(1) * (-t23 * pkin(1) + t26) - g(2) * (t25 * pkin(1) + t27); 0, 0, 0, 0, t8, t9, 0, 0, 0, 0, 0, t6, -t5, 0, 0, 0, 0, 0, t4, -t3, t4, -t3, -t9, -g(1) * t26 - g(2) * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t24 + t9 * t22, g(3) * t22 + t9 * t24, 0, 0, 0, 0, 0, t1, t2, t1, t2, 0, -g(3) * t28 - t9 * (-t22 * pkin(3) - pkin(4) * t14); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8;];
taug_reg = t7;
