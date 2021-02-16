% Calculate minimal parameter regressor of gravitation load for
% S5RRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% taug_reg [5x28]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:50
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:48:58
% EndTime: 2021-01-15 22:48:58
% DurationCPUTime: 0.11s
% Computational Cost: add. (184->27), mult. (149->41), div. (0->0), fcn. (141->10), ass. (0->25)
t22 = qJ(2) + qJ(3);
t19 = cos(t22);
t25 = cos(qJ(2));
t27 = t25 * pkin(2) + pkin(3) * t19;
t17 = pkin(9) + t22;
t24 = sin(qJ(1));
t26 = cos(qJ(1));
t10 = g(1) * t26 + g(2) * t24;
t9 = g(1) * t24 - g(2) * t26;
t18 = sin(t22);
t5 = -g(3) * t19 + t10 * t18;
t23 = sin(qJ(2));
t21 = -qJ(4) - pkin(7) - pkin(6);
t16 = qJ(5) + t17;
t14 = cos(t17);
t13 = sin(t17);
t12 = cos(t16);
t11 = sin(t16);
t7 = pkin(1) + t27;
t6 = g(3) * t18 + t10 * t19;
t4 = g(3) * t13 + t10 * t14;
t3 = -g(3) * t14 + t10 * t13;
t2 = g(3) * t11 + t10 * t12;
t1 = -g(3) * t12 + t10 * t11;
t8 = [0, t9, t10, 0, 0, 0, 0, 0, t9 * t25, -t9 * t23, 0, 0, 0, 0, 0, t9 * t19, -t9 * t18, t9 * t14, -t9 * t13, -t10, -g(1) * (-t26 * t21 - t24 * t7) - g(2) * (-t24 * t21 + t26 * t7), 0, 0, 0, 0, 0, t9 * t12, -t9 * t11; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t25 + t10 * t23, g(3) * t23 + t10 * t25, 0, 0, 0, 0, 0, t5, t6, t3, t4, 0, -g(3) * t27 - t10 * (-t23 * pkin(2) - pkin(3) * t18), 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, t3, t4, 0, t5 * pkin(3), 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t8;
