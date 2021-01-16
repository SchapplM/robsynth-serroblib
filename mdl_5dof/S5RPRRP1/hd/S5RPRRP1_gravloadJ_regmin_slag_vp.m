% Calculate minimal parameter regressor of gravitation load for
% S5RPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% taug_reg [5x24]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:27
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:27:19
% EndTime: 2021-01-15 12:27:20
% DurationCPUTime: 0.11s
% Computational Cost: add. (97->28), mult. (129->35), div. (0->0), fcn. (119->6), ass. (0->17)
t13 = sin(qJ(3));
t12 = qJ(3) + qJ(4);
t8 = sin(t12);
t24 = -t13 * pkin(3) - pkin(4) * t8;
t14 = sin(qJ(1));
t16 = cos(qJ(1));
t23 = -g(1) * t14 + g(2) * t16;
t18 = qJ(2) - t24;
t6 = g(1) * t16 + g(2) * t14;
t9 = cos(t12);
t2 = g(3) * t8 + t23 * t9;
t15 = cos(qJ(3));
t10 = qJ(5) + pkin(1) + pkin(6) + pkin(7);
t4 = t6 * t9;
t3 = t6 * t8;
t1 = g(3) * t9 - t23 * t8;
t5 = [0, -t23, t6, t23, -t6, -g(1) * (-t14 * pkin(1) + t16 * qJ(2)) - g(2) * (t16 * pkin(1) + t14 * qJ(2)), 0, 0, 0, 0, 0, -t6 * t13, -t6 * t15, 0, 0, 0, 0, 0, -t3, -t4, -t3, -t4, -t23, -g(1) * (-t10 * t14 + t18 * t16) - g(2) * (t10 * t16 + t18 * t14); 0, 0, 0, 0, 0, t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t13 + t15 * t23, g(3) * t15 - t13 * t23, 0, 0, 0, 0, 0, t2, t1, t2, t1, 0, -g(3) * t24 + t23 * (pkin(3) * t15 + pkin(4) * t9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, t2, t1, 0, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6;];
taug_reg = t5;
