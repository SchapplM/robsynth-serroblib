% Calculate minimal parameter regressor of gravitation load for
% S4RRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% taug_reg [4x24]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRRR4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR4_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR4_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:26:14
% EndTime: 2019-12-31 17:26:15
% DurationCPUTime: 0.12s
% Computational Cost: add. (78->21), mult. (118->36), div. (0->0), fcn. (126->8), ass. (0->25)
t11 = qJ(2) + qJ(3);
t10 = cos(t11);
t14 = sin(qJ(1));
t17 = cos(qJ(1));
t19 = g(1) * t17 + g(2) * t14;
t9 = sin(t11);
t3 = -g(3) * t10 + t19 * t9;
t25 = g(3) * t9;
t12 = sin(qJ(4));
t23 = t14 * t12;
t15 = cos(qJ(4));
t22 = t14 * t15;
t21 = t17 * t12;
t20 = t17 * t15;
t18 = g(1) * t14 - g(2) * t17;
t16 = cos(qJ(2));
t13 = sin(qJ(2));
t8 = t10 * t20 + t23;
t7 = -t10 * t21 + t22;
t6 = -t10 * t22 + t21;
t5 = t10 * t23 + t20;
t4 = t19 * t10 + t25;
t2 = t3 * t15;
t1 = t3 * t12;
t24 = [0, t18, t19, 0, 0, 0, 0, 0, t18 * t16, -t18 * t13, 0, 0, 0, 0, 0, t18 * t10, -t18 * t9, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t8, -g(1) * t5 - g(2) * t7; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t16 + t19 * t13, g(3) * t13 + t19 * t16, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 + g(2) * t5 + t12 * t25, g(1) * t8 - g(2) * t6 + t15 * t25;];
taug_reg = t24;
