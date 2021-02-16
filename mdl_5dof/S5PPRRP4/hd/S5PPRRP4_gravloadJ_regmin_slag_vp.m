% Calculate minimal parameter regressor of gravitation load for
% S5PPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
% 
% Output:
% taug_reg [5x16]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:56
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PPRRP4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:56:28
% EndTime: 2021-01-15 14:56:29
% DurationCPUTime: 0.08s
% Computational Cost: add. (60->19), mult. (124->24), div. (0->0), fcn. (149->6), ass. (0->18)
t20 = cos(qJ(3));
t19 = sin(qJ(3));
t18 = cos(pkin(7));
t17 = sin(pkin(7));
t5 = -t17 * t19 - t18 * t20;
t6 = -t17 * t20 + t18 * t19;
t16 = g(1) * t6 - g(2) * t5;
t15 = g(1) * t5 + g(2) * t6;
t13 = sin(qJ(4));
t14 = cos(qJ(4));
t1 = g(3) * t14 - t15 * t13;
t12 = qJ(5) + pkin(6);
t11 = t14 * pkin(4) + pkin(3);
t7 = -g(1) * t17 + g(2) * t18;
t4 = t16 * t14;
t3 = t16 * t13;
t2 = -g(3) * t13 - t15 * t14;
t8 = [-g(3), -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, t16, -t15, 0, 0, 0, 0, 0, t4, -t3, t4, -t3, t15, -g(1) * (-t6 * t11 - t5 * t12) - g(2) * (t5 * t11 - t6 * t12); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16;];
taug_reg = t8;
