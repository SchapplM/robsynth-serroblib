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
% taug_reg [5x20]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 16:19:40
% EndTime: 2019-12-05 16:19:41
% DurationCPUTime: 0.08s
% Computational Cost: add. (97->19), mult. (70->24), div. (0->0), fcn. (65->6), ass. (0->16)
t11 = pkin(8) + qJ(2);
t8 = sin(t11);
t9 = cos(t11);
t4 = g(1) * t9 + g(2) * t8;
t3 = g(1) * t8 - g(2) * t9;
t13 = sin(qJ(3));
t14 = cos(qJ(3));
t15 = -g(3) * t14 + t4 * t13;
t12 = -qJ(4) - pkin(6);
t10 = qJ(3) + pkin(9) + qJ(5);
t7 = t14 * pkin(3) + pkin(2);
t6 = cos(t10);
t5 = sin(t10);
t2 = g(3) * t5 + t4 * t6;
t1 = -g(3) * t6 + t4 * t5;
t16 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, t3, t4, 0, 0, 0, 0, 0, t3 * t14, -t3 * t13, -t4, -g(1) * (-t9 * t12 - t8 * t7) - g(2) * (-t8 * t12 + t9 * t7), 0, 0, 0, 0, 0, t3 * t6, -t3 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, g(3) * t13 + t4 * t14, 0, t15 * pkin(3), 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t16;
