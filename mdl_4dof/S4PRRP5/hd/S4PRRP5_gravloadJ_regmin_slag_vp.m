% Calculate minimal parameter regressor of gravitation load for
% S4PRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% taug_reg [4x15]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-14 22:36
% Revision: beb2ba9bd8c5bd556f42a244985f3dab86917626 (2021-01-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4PRRP5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-14 22:36:08
% EndTime: 2021-01-14 22:36:09
% DurationCPUTime: 0.11s
% Computational Cost: add. (49->22), mult. (123->32), div. (0->0), fcn. (125->6), ass. (0->19)
t8 = sin(pkin(6));
t9 = cos(pkin(6));
t23 = -g(1) * t9 - g(2) * t8;
t12 = sin(qJ(2));
t14 = cos(qJ(2));
t5 = -g(3) * t14 - t23 * t12;
t20 = g(3) * t12;
t11 = sin(qJ(3));
t18 = t11 * t14;
t13 = cos(qJ(3));
t17 = t13 * t14;
t1 = -g(1) * (t8 * t13 - t9 * t18) - g(2) * (-t9 * t13 - t8 * t18) + t11 * t20;
t10 = qJ(4) + pkin(5);
t7 = t13 * pkin(3) + pkin(2);
t6 = -t23 * t14 + t20;
t4 = t5 * t13;
t3 = t5 * t11;
t2 = -g(1) * (-t8 * t11 - t9 * t17) - g(2) * (t9 * t11 - t8 * t17) + t13 * t20;
t15 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, t5, t6, 0, 0, 0, 0, 0, t4, -t3, t4, -t3, -t6, -g(3) * (t12 * t10 + t14 * t7) + t23 * (t10 * t14 - t12 * t7); 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, t2, 0, t1 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5;];
taug_reg = t15;
