% Calculate minimal parameter regressor of gravitation load for
% S5RRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% taug_reg [5x24]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPR4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:11:39
% EndTime: 2019-12-31 21:11:40
% DurationCPUTime: 0.16s
% Computational Cost: add. (169->34), mult. (208->48), div. (0->0), fcn. (221->8), ass. (0->32)
t26 = sin(qJ(3));
t29 = cos(qJ(3));
t36 = t29 * pkin(3) + t26 * qJ(4);
t46 = -pkin(2) - t36;
t25 = sin(qJ(5));
t28 = cos(qJ(5));
t34 = t29 * t25 - t26 * t28;
t24 = qJ(1) + qJ(2);
t23 = cos(t24);
t22 = sin(t24);
t45 = g(1) * t22;
t11 = -g(2) * t23 + t45;
t12 = g(1) * t23 + g(2) * t22;
t37 = t22 * pkin(7) - t46 * t23;
t13 = t26 * t25 + t29 * t28;
t5 = t34 * t22;
t7 = t34 * t23;
t33 = g(1) * t7 + g(2) * t5 + g(3) * t13;
t6 = t13 * t22;
t8 = t13 * t23;
t32 = g(1) * t8 + g(2) * t6 - g(3) * t34;
t31 = t46 * t45;
t30 = cos(qJ(1));
t27 = sin(qJ(1));
t20 = t23 * pkin(7);
t10 = t11 * t29;
t9 = t11 * t26;
t4 = g(3) * t26 + t12 * t29;
t3 = -g(3) * t29 + t12 * t26;
t2 = g(1) * t6 - g(2) * t8;
t1 = -g(1) * t5 + g(2) * t7;
t14 = [0, g(1) * t27 - g(2) * t30, g(1) * t30 + g(2) * t27, 0, t11, t12, 0, 0, 0, 0, 0, t10, -t9, t10, -t12, t9, -g(1) * (-t27 * pkin(1) + t20) - g(2) * (t30 * pkin(1) + t37) - t31, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, t11, t12, 0, 0, 0, 0, 0, t10, -t9, t10, -t12, t9, -g(1) * t20 - g(2) * t37 - t31, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, t3, 0, -t4, -g(3) * t36 + t12 * (pkin(3) * t26 - qJ(4) * t29), 0, 0, 0, 0, 0, -t33, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, t32;];
taug_reg = t14;
