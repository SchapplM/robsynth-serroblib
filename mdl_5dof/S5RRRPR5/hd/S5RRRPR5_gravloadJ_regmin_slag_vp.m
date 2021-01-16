% Calculate minimal parameter regressor of gravitation load for
% S5RRRPR5
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
% Datum: 2021-01-15 23:12
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 23:10:44
% EndTime: 2021-01-15 23:10:44
% DurationCPUTime: 0.13s
% Computational Cost: add. (170->35), mult. (179->55), div. (0->0), fcn. (179->10), ass. (0->33)
t26 = sin(qJ(1));
t29 = cos(qJ(1));
t14 = g(1) * t29 + g(2) * t26;
t23 = qJ(2) + qJ(3);
t18 = pkin(9) + t23;
t15 = sin(t18);
t16 = cos(t18);
t3 = -g(3) * t16 + t14 * t15;
t36 = g(3) * t15;
t24 = sin(qJ(5));
t34 = t26 * t24;
t27 = cos(qJ(5));
t33 = t26 * t27;
t32 = t29 * t24;
t31 = t29 * t27;
t20 = cos(t23);
t28 = cos(qJ(2));
t30 = t28 * pkin(2) + pkin(3) * t20;
t13 = g(1) * t26 - g(2) * t29;
t19 = sin(t23);
t5 = -g(3) * t20 + t14 * t19;
t25 = sin(qJ(2));
t22 = -qJ(4) - pkin(7) - pkin(6);
t11 = pkin(1) + t30;
t10 = t16 * t31 + t34;
t9 = -t16 * t32 + t33;
t8 = -t16 * t33 + t32;
t7 = t16 * t34 + t31;
t6 = g(3) * t19 + t14 * t20;
t4 = t14 * t16 + t36;
t2 = t3 * t27;
t1 = t3 * t24;
t12 = [0, t13, t14, 0, 0, 0, 0, 0, t13 * t28, -t13 * t25, 0, 0, 0, 0, 0, t13 * t20, -t13 * t19, t13 * t16, -t13 * t15, -t14, -g(1) * (-t26 * t11 - t29 * t22) - g(2) * (t29 * t11 - t26 * t22), 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t10, -g(1) * t7 - g(2) * t9; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t28 + t14 * t25, g(3) * t25 + t14 * t28, 0, 0, 0, 0, 0, t5, t6, t3, t4, 0, -g(3) * t30 - t14 * (-t25 * pkin(2) - pkin(3) * t19), 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, t3, t4, 0, t5 * pkin(3), 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 + g(2) * t7 + t24 * t36, g(1) * t10 - g(2) * t8 + t27 * t36;];
taug_reg = t12;
