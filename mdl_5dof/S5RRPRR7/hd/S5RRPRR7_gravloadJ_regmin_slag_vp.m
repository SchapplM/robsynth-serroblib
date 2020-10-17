% Calculate minimal parameter regressor of gravitation load for
% S5RRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% taug_reg [5x23]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:15:40
% EndTime: 2019-12-31 20:15:40
% DurationCPUTime: 0.09s
% Computational Cost: add. (113->29), mult. (98->30), div. (0->0), fcn. (94->8), ass. (0->21)
t18 = qJ(1) + qJ(2);
t14 = sin(t18);
t16 = cos(t18);
t24 = t16 * pkin(2) + t14 * qJ(3);
t23 = -t14 * pkin(2) + t16 * qJ(3);
t8 = g(1) * t16 + g(2) * t14;
t7 = g(1) * t14 - g(2) * t16;
t22 = cos(qJ(1));
t21 = cos(qJ(4));
t20 = sin(qJ(1));
t19 = sin(qJ(4));
t17 = qJ(4) + qJ(5);
t15 = cos(t17);
t13 = sin(t17);
t6 = t8 * t21;
t5 = t8 * t19;
t4 = t8 * t15;
t3 = t8 * t13;
t2 = g(3) * t13 - t7 * t15;
t1 = g(3) * t15 + t7 * t13;
t9 = [0, g(1) * t20 - g(2) * t22, g(1) * t22 + g(2) * t20, 0, t7, t8, -t7, -t8, -g(1) * (-t20 * pkin(1) + t23) - g(2) * (t22 * pkin(1) + t24), 0, 0, 0, 0, 0, -t5, -t6, 0, 0, 0, 0, 0, -t3, -t4; 0, 0, 0, 0, t7, t8, -t7, -t8, -g(1) * t23 - g(2) * t24, 0, 0, 0, 0, 0, -t5, -t6, 0, 0, 0, 0, 0, -t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t19 - t7 * t21, g(3) * t21 + t7 * t19, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1;];
taug_reg = t9;
