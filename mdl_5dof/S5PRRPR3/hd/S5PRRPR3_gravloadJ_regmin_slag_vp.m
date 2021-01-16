% Calculate minimal parameter regressor of gravitation load for
% S5PRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% taug_reg [5x22]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:42
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRPR3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:42:09
% EndTime: 2021-01-15 15:42:10
% DurationCPUTime: 0.11s
% Computational Cost: add. (121->23), mult. (88->30), div. (0->0), fcn. (83->8), ass. (0->19)
t14 = qJ(3) + pkin(9);
t13 = pkin(8) + qJ(2);
t10 = cos(t13);
t8 = sin(t13);
t3 = g(1) * t8 - g(2) * t10;
t4 = g(1) * t10 + g(2) * t8;
t16 = sin(qJ(3));
t17 = cos(qJ(3));
t18 = -g(3) * t17 + t4 * t16;
t15 = -qJ(4) - pkin(6);
t12 = qJ(5) + t14;
t11 = cos(t14);
t9 = sin(t14);
t7 = t17 * pkin(3) + pkin(2);
t6 = cos(t12);
t5 = sin(t12);
t2 = g(3) * t5 + t4 * t6;
t1 = -g(3) * t6 + t4 * t5;
t19 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, t3, t4, 0, 0, 0, 0, 0, t3 * t17, -t3 * t16, t3 * t11, -t3 * t9, -t4, -g(1) * (-t10 * t15 - t8 * t7) - g(2) * (t10 * t7 - t8 * t15), 0, 0, 0, 0, 0, t3 * t6, -t3 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, g(3) * t16 + t4 * t17, -g(3) * t11 + t4 * t9, g(3) * t9 + t4 * t11, 0, t18 * pkin(3), 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t19;
