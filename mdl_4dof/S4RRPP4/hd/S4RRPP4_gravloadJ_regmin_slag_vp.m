% Calculate minimal parameter regressor of gravitation load for
% S4RRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% taug_reg [4x18]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRPP4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_gravloadJ_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:59:18
% EndTime: 2019-12-31 16:59:18
% DurationCPUTime: 0.12s
% Computational Cost: add. (58->32), mult. (139->41), div. (0->0), fcn. (129->4), ass. (0->24)
t17 = sin(qJ(2));
t18 = sin(qJ(1));
t20 = cos(qJ(1));
t6 = g(1) * t20 + g(2) * t18;
t32 = t6 * t17;
t11 = t17 * qJ(3);
t19 = cos(qJ(2));
t24 = t19 * pkin(2) + t11;
t30 = pkin(2) * t17;
t29 = g(1) * t18;
t26 = t19 * pkin(3);
t25 = t19 * t20;
t23 = qJ(3) * t19;
t22 = pkin(2) * t25 + t18 * pkin(5) + (pkin(1) + t11) * t20;
t5 = -g(2) * t20 + t29;
t21 = -pkin(1) - t24;
t14 = t20 * pkin(5);
t9 = t20 * t23;
t7 = t18 * t23;
t4 = t5 * t19;
t3 = t5 * t17;
t2 = g(3) * t17 + t6 * t19;
t1 = -g(3) * t19 + t32;
t8 = [0, t5, t6, 0, 0, 0, 0, 0, t4, -t3, t4, -t6, t3, -g(1) * t14 - g(2) * t22 - t21 * t29, t4, t3, t6, -g(1) * (-t20 * qJ(4) + t14) - g(2) * (pkin(3) * t25 + t22) + (-g(1) * (t21 - t26) + g(2) * qJ(4)) * t18; 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, 0, -t2, -g(1) * (-t20 * t30 + t9) - g(2) * (-t18 * t30 + t7) - g(3) * t24, t1, -t2, 0, -g(1) * t9 - g(2) * t7 - g(3) * (t24 + t26) + (pkin(2) + pkin(3)) * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5;];
taug_reg = t8;
