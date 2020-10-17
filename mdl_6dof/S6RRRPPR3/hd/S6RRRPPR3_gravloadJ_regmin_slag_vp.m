% Calculate minimal parameter regressor of gravitation load for
% S6RRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPPR3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 04:32:42
% EndTime: 2019-05-07 04:32:43
% DurationCPUTime: 0.28s
% Computational Cost: add. (241->64), mult. (300->75), div. (0->0), fcn. (292->8), ass. (0->44)
t31 = cos(qJ(2));
t26 = qJ(2) + qJ(3);
t23 = sin(t26);
t24 = cos(t26);
t41 = t24 * pkin(3) + t23 * qJ(4);
t38 = t31 * pkin(2) + t41;
t55 = pkin(1) + t38;
t32 = cos(qJ(1));
t54 = g(2) * t32;
t29 = sin(qJ(1));
t49 = g(1) * t32;
t12 = g(2) * t29 + t49;
t53 = t12 * t23;
t4 = g(3) * t23 + t12 * t24;
t28 = sin(qJ(2));
t51 = pkin(2) * t28;
t50 = pkin(3) * t23;
t46 = g(3) * t24;
t19 = t24 * pkin(4);
t27 = sin(qJ(6));
t45 = t29 * t27;
t30 = cos(qJ(6));
t44 = t29 * t30;
t43 = t32 * t27;
t42 = t32 * t30;
t40 = qJ(4) * t24;
t33 = -pkin(8) - pkin(7);
t39 = qJ(5) + t33;
t37 = t55 * t54;
t36 = -t50 - t51;
t11 = g(1) * t29 - t54;
t34 = (pkin(3) + pkin(4)) * t53;
t15 = t32 * t40;
t13 = t29 * t40;
t10 = t23 * t42 - t45;
t9 = -t23 * t43 - t44;
t8 = -t23 * t44 - t43;
t7 = t23 * t45 - t42;
t6 = t11 * t24;
t5 = t11 * t23;
t3 = -t46 + t53;
t2 = t4 * t30;
t1 = t4 * t27;
t14 = [0, t11, t12, 0, 0, 0, 0, 0, t11 * t31, -t11 * t28, 0, 0, 0, 0, 0, t6, -t5, t6, -t12, t5, t33 * t49 - t37 + (g(1) * t55 + g(2) * t33) * t29, t5, -t6, t12, -t37 + (g(1) * t39 - g(2) * t19) * t32 + (-g(1) * (-t55 - t19) + g(2) * t39) * t29, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t10, -g(1) * t7 - g(2) * t9; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t31 + t12 * t28, g(3) * t28 + t12 * t31, 0, 0, 0, 0, 0, t3, t4, t3, 0, -t4, -g(1) * (t36 * t32 + t15) - g(2) * (t36 * t29 + t13) - g(3) * t38, -t4, -t3, 0, -g(1) * (-t32 * t51 + t15) - g(2) * (-t29 * t51 + t13) - g(3) * (t19 + t38) + t34, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, t3, 0, -t4, -g(1) * (-t32 * t50 + t15) - g(2) * (-t29 * t50 + t13) - g(3) * t41, -t4, -t3, 0, -g(1) * t15 - g(2) * t13 - g(3) * (t19 + t41) + t34, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 + g(2) * t7 - t27 * t46, g(1) * t10 - g(2) * t8 - t30 * t46;];
taug_reg  = t14;
