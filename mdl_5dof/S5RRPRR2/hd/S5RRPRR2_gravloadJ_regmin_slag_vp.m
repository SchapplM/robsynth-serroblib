% Calculate minimal parameter regressor of gravitation load for
% S5RRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% taug_reg [5x28]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:24
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:23:35
% EndTime: 2021-01-15 21:23:36
% DurationCPUTime: 0.10s
% Computational Cost: add. (157->24), mult. (126->36), div. (0->0), fcn. (121->10), ass. (0->23)
t16 = qJ(2) + pkin(9);
t15 = qJ(4) + t16;
t19 = sin(qJ(1));
t21 = cos(qJ(1));
t6 = g(1) * t21 + g(2) * t19;
t5 = g(1) * t19 - g(2) * t21;
t18 = sin(qJ(2));
t20 = cos(qJ(2));
t22 = -g(3) * t20 + t6 * t18;
t17 = -qJ(3) - pkin(6);
t14 = cos(t16);
t13 = sin(t16);
t12 = t20 * pkin(2) + pkin(1);
t11 = qJ(5) + t15;
t10 = cos(t15);
t9 = sin(t15);
t8 = cos(t11);
t7 = sin(t11);
t4 = g(3) * t9 + t6 * t10;
t3 = -g(3) * t10 + t6 * t9;
t2 = g(3) * t7 + t6 * t8;
t1 = -g(3) * t8 + t6 * t7;
t23 = [0, t5, t6, 0, 0, 0, 0, 0, t5 * t20, -t5 * t18, t5 * t14, -t5 * t13, -t6, -g(1) * (-t19 * t12 - t17 * t21) - g(2) * (t21 * t12 - t19 * t17), 0, 0, 0, 0, 0, t5 * t10, -t5 * t9, 0, 0, 0, 0, 0, t5 * t8, -t5 * t7; 0, 0, 0, 0, 0, 0, 0, 0, t22, g(3) * t18 + t6 * t20, -g(3) * t14 + t6 * t13, g(3) * t13 + t6 * t14, 0, t22 * pkin(2), 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t23;
