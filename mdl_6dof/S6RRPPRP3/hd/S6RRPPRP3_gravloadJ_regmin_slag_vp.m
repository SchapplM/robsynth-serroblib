% Calculate minimal parameter regressor of gravitation load for
% S6RRPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
% 
% Output:
% taug_reg [6x27]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPPRP3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:16:35
% EndTime: 2019-05-06 09:16:36
% DurationCPUTime: 0.28s
% Computational Cost: add. (129->67), mult. (292->84), div. (0->0), fcn. (282->6), ass. (0->50)
t26 = sin(qJ(2));
t17 = t26 * qJ(3);
t61 = pkin(1) + t17;
t30 = cos(qJ(1));
t27 = sin(qJ(1));
t51 = g(2) * t27;
t10 = g(1) * t30 + t51;
t60 = t10 * t26;
t25 = sin(qJ(5));
t28 = cos(qJ(5));
t42 = t30 * t28;
t47 = t27 * t25;
t3 = t26 * t47 - t42;
t29 = cos(qJ(2));
t49 = g(3) * t29;
t43 = t30 * t25;
t46 = t27 * t28;
t5 = -t26 * t43 - t46;
t59 = -g(1) * t5 + g(2) * t3 - t25 * t49;
t2 = g(3) * t26 + t10 * t29;
t56 = pkin(2) + pkin(3);
t55 = pkin(2) * t26;
t21 = t30 * pkin(7);
t54 = g(1) * t21;
t53 = g(1) * t27;
t20 = t29 * pkin(2);
t19 = t29 * pkin(3);
t16 = t28 * pkin(5) + pkin(4);
t48 = t26 * t16;
t24 = -qJ(6) - pkin(8);
t45 = t29 * t24;
t44 = t29 * t30;
t41 = t20 + t17;
t40 = qJ(3) * t29;
t39 = -t24 + t56;
t37 = t19 + t41;
t35 = pkin(5) * t25 + qJ(4);
t34 = pkin(2) * t44 + t27 * pkin(7) + t61 * t30;
t33 = g(1) * t39;
t9 = -g(2) * t30 + t53;
t32 = -t61 - t20;
t31 = g(2) * (pkin(3) * t44 + t34);
t13 = t30 * t40;
t11 = t27 * t40;
t8 = t9 * t29;
t7 = t9 * t26;
t6 = t26 * t42 - t47;
t4 = -t26 * t46 - t43;
t1 = -t49 + t60;
t12 = [0, t9, t10, 0, 0, 0, 0, 0, t8, -t7, t8, -t10, t7, -g(2) * t34 - t32 * t53 - t54, t7, -t8, t10, -g(1) * (-t30 * qJ(4) + t21) - t31 + (-g(1) * (t32 - t19) + g(2) * qJ(4)) * t27, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5, t8, -t54 - t31 + (g(1) * t35 - g(2) * (-t45 + t48)) * t30 + (-g(1) * (-t61 - t48) + g(2) * t35 + t29 * t33) * t27; 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, 0, -t2, -g(1) * (-t30 * t55 + t13) - g(2) * (-t27 * t55 + t11) - g(3) * t41, -t2, -t1, 0, -g(1) * t13 - g(2) * t11 - g(3) * t37 + t56 * t60, 0, 0, 0, 0, 0, -t2 * t28, t2 * t25, t1, -g(1) * (t16 * t44 + t13) - g(2) * (t27 * t29 * t16 + t11) - g(3) * (t37 - t45) + (-g(3) * t16 + t30 * t33 + t39 * t51) * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, g(1) * t6 - g(2) * t4 - t28 * t49, 0, t59 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2;];
taug_reg  = t12;
