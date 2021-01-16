% Calculate minimal parameter regressor of gravitation load for
% S4RRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% taug_reg [4x17]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:04
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRRP2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:04:29
% EndTime: 2021-01-15 11:04:29
% DurationCPUTime: 0.08s
% Computational Cost: add. (91->22), mult. (96->26), div. (0->0), fcn. (89->6), ass. (0->18)
t11 = qJ(1) + qJ(2);
t10 = cos(t11);
t12 = -qJ(4) - pkin(6);
t15 = cos(qJ(3));
t8 = t15 * pkin(3) + pkin(2);
t9 = sin(t11);
t18 = t10 * t8 - t9 * t12;
t5 = g(1) * t9 - g(2) * t10;
t6 = g(1) * t10 + g(2) * t9;
t17 = -t10 * t12 - t9 * t8;
t13 = sin(qJ(3));
t1 = -g(3) * t15 + t6 * t13;
t16 = cos(qJ(1));
t14 = sin(qJ(1));
t4 = t5 * t15;
t3 = t5 * t13;
t2 = g(3) * t13 + t6 * t15;
t7 = [0, g(1) * t14 - g(2) * t16, g(1) * t16 + g(2) * t14, 0, t5, t6, 0, 0, 0, 0, 0, t4, -t3, t4, -t3, -t6, -g(1) * (-t14 * pkin(1) + t17) - g(2) * (t16 * pkin(1) + t18); 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t4, -t3, t4, -t3, -t6, -g(1) * t17 - g(2) * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, t2, 0, t1 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5;];
taug_reg = t7;
