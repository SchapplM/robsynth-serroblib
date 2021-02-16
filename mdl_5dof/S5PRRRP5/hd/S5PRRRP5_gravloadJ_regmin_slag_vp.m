% Calculate minimal parameter regressor of gravitation load for
% S5PRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% taug_reg [5x22]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:34
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRRP5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:33:42
% EndTime: 2021-01-15 16:33:43
% DurationCPUTime: 0.15s
% Computational Cost: add. (157->38), mult. (222->60), div. (0->0), fcn. (233->8), ass. (0->26)
t16 = sin(pkin(8));
t17 = cos(pkin(8));
t24 = g(1) * t17 + g(2) * t16;
t19 = sin(qJ(2));
t21 = cos(qJ(2));
t5 = -g(3) * t21 + t24 * t19;
t30 = g(3) * t19;
t28 = t16 * t21;
t27 = t17 * t21;
t18 = sin(qJ(3));
t26 = t18 * t21;
t20 = cos(qJ(3));
t25 = t20 * t21;
t15 = qJ(3) + qJ(4);
t12 = cos(t15);
t9 = t20 * pkin(3) + pkin(4) * t12;
t11 = sin(t15);
t1 = -g(1) * (-t11 * t27 + t16 * t12) - g(2) * (-t11 * t28 - t17 * t12) + t11 * t30;
t14 = -qJ(5) - pkin(7) - pkin(6);
t8 = -t18 * pkin(3) - pkin(4) * t11;
t7 = pkin(2) + t9;
t6 = t24 * t21 + t30;
t4 = t5 * t12;
t3 = t5 * t11;
t2 = -g(1) * (-t16 * t11 - t12 * t27) - g(2) * (t17 * t11 - t12 * t28) + t12 * t30;
t10 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, t5, t6, 0, 0, 0, 0, 0, t5 * t20, -t5 * t18, 0, 0, 0, 0, 0, t4, -t3, t4, -t3, -t6, -g(3) * (-t19 * t14 + t21 * t7) + t24 * (t14 * t21 + t19 * t7); 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t16 * t20 - t17 * t26) - g(2) * (-t16 * t26 - t17 * t20) + t18 * t30, -g(1) * (-t16 * t18 - t17 * t25) - g(2) * (-t16 * t25 + t17 * t18) + t20 * t30, 0, 0, 0, 0, 0, t1, t2, t1, t2, 0, -g(1) * (t16 * t9 + t8 * t27) - g(2) * (-t17 * t9 + t8 * t28) - t8 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5;];
taug_reg = t10;
