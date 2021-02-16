% Calculate minimal parameter regressor of gravitation load for
% S5PRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% taug_reg [5x18]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:05
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRPRP2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:04:51
% EndTime: 2021-01-15 15:04:52
% DurationCPUTime: 0.12s
% Computational Cost: add. (126->32), mult. (136->47), div. (0->0), fcn. (145->6), ass. (0->25)
t19 = sin(pkin(8));
t20 = cos(pkin(8));
t23 = cos(qJ(4));
t34 = -t19 * (-qJ(5) - pkin(6)) + (t23 * pkin(4) + pkin(3)) * t20;
t22 = sin(qJ(4));
t31 = g(3) * t19;
t18 = pkin(7) + qJ(2);
t16 = sin(t18);
t17 = cos(t18);
t27 = t20 * t22;
t5 = t16 * t27 + t17 * t23;
t7 = t16 * t23 - t17 * t27;
t1 = -g(1) * t7 + g(2) * t5 + t22 * t31;
t29 = t17 * t22;
t26 = t20 * t23;
t25 = t17 * pkin(2) + t16 * qJ(3);
t10 = g(1) * t17 + g(2) * t16;
t9 = g(1) * t16 - g(2) * t17;
t12 = t17 * qJ(3);
t8 = t16 * t22 + t17 * t26;
t6 = -t16 * t26 + t29;
t4 = -g(1) * t6 - g(2) * t8;
t3 = -g(1) * t5 - g(2) * t7;
t2 = g(1) * t8 - g(2) * t6 + t23 * t31;
t11 = [-g(3), 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, t9, t10, t9 * t20, -t10, -g(1) * (-t16 * pkin(2) + t12) - g(2) * t25, 0, 0, 0, 0, 0, t4, t3, t4, t3, t9 * t19, -g(1) * (pkin(4) * t29 + t12) - g(2) * (t34 * t17 + t25) + (-g(1) * (-pkin(2) - t34) - g(2) * pkin(4) * t22) * t16; 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t20 - t10 * t19;];
taug_reg = t11;
