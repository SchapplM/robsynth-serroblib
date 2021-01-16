% Calculate minimal parameter regressor of gravitation load for
% S5RPRRP12
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
% Datum: 2021-01-15 19:26
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRP12_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP12_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:25:53
% EndTime: 2021-01-15 19:25:54
% DurationCPUTime: 0.14s
% Computational Cost: add. (86->36), mult. (197->51), div. (0->0), fcn. (203->6), ass. (0->32)
t24 = cos(qJ(4));
t18 = t24 * pkin(4) + pkin(3);
t20 = -qJ(5) - pkin(7);
t22 = sin(qJ(3));
t25 = cos(qJ(3));
t43 = -t18 * t22 - t20 * t25;
t42 = qJ(2) - t43;
t23 = sin(qJ(1));
t26 = cos(qJ(1));
t41 = -g(1) * t23 + g(2) * t26;
t21 = sin(qJ(4));
t31 = t26 * t21;
t32 = t23 * t24;
t11 = t22 * t31 + t32;
t35 = g(3) * t25;
t30 = t26 * t24;
t33 = t23 * t21;
t9 = -t22 * t33 + t30;
t1 = -g(1) * t9 - g(2) * t11 + t21 * t35;
t8 = -g(3) * t22 - t25 * t41;
t15 = g(1) * t26 + g(2) * t23;
t16 = t21 * pkin(4) + pkin(1) + pkin(6);
t13 = t15 * t25;
t12 = t22 * t30 - t33;
t10 = t22 * t32 + t31;
t7 = -t22 * t41 + t35;
t6 = t8 * t24;
t5 = t8 * t21;
t4 = -g(1) * t12 - g(2) * t10;
t3 = g(1) * t11 - g(2) * t9;
t2 = g(1) * t10 - g(2) * t12 + t24 * t35;
t14 = [0, -t41, t15, t41, -t15, -g(1) * (-t23 * pkin(1) + t26 * qJ(2)) - g(2) * (t26 * pkin(1) + t23 * qJ(2)), 0, 0, 0, 0, 0, -t15 * t22, -t13, 0, 0, 0, 0, 0, t4, t3, t4, t3, t13, -g(1) * (-t16 * t23 + t42 * t26) - g(2) * (t16 * t26 + t42 * t23); 0, 0, 0, 0, 0, t41, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, t7, 0, 0, 0, 0, 0, -t6, t5, -t6, t5, -t7, -g(3) * t43 + t41 * (t25 * t18 - t20 * t22); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8;];
taug_reg = t14;
