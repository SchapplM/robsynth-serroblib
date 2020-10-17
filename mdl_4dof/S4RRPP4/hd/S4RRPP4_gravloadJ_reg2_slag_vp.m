% Calculate inertial parameters regressor of gravitation load for
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
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRPP4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_gravloadJ_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:59:18
% EndTime: 2019-12-31 16:59:19
% DurationCPUTime: 0.13s
% Computational Cost: add. (62->35), mult. (147->45), div. (0->0), fcn. (135->4), ass. (0->25)
t17 = sin(qJ(2));
t18 = sin(qJ(1));
t20 = cos(qJ(1));
t6 = g(1) * t20 + g(2) * t18;
t33 = t6 * t17;
t11 = t17 * qJ(3);
t19 = cos(qJ(2));
t25 = t19 * pkin(2) + t11;
t31 = pkin(2) * t17;
t30 = g(1) * t18;
t27 = t19 * pkin(3);
t26 = t19 * t20;
t24 = t20 * pkin(1) + t18 * pkin(5);
t23 = qJ(3) * t19;
t22 = pkin(2) * t26 + t20 * t11 + t24;
t5 = -g(2) * t20 + t30;
t21 = -pkin(1) - t25;
t14 = t20 * pkin(5);
t9 = t20 * t23;
t7 = t18 * t23;
t4 = t5 * t19;
t3 = t5 * t17;
t2 = g(3) * t17 + t6 * t19;
t1 = -g(3) * t19 + t33;
t8 = [0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (-t18 * pkin(1) + t14) - g(2) * t24, 0, 0, 0, 0, 0, 0, t4, -t6, t3, -g(1) * t14 - g(2) * t22 - t21 * t30, 0, 0, 0, 0, 0, 0, t4, t3, t6, -g(1) * (-t20 * qJ(4) + t14) - g(2) * (pkin(3) * t26 + t22) + (-g(1) * (t21 - t27) + g(2) * qJ(4)) * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t2, -g(1) * (-t20 * t31 + t9) - g(2) * (-t18 * t31 + t7) - g(3) * t25, 0, 0, 0, 0, 0, 0, t1, -t2, 0, -g(1) * t9 - g(2) * t7 - g(3) * (t25 + t27) + (pkin(2) + pkin(3)) * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5;];
taug_reg = t8;
