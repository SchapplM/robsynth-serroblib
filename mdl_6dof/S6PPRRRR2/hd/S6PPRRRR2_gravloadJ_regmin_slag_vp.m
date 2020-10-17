% Calculate minimal parameter regressor of gravitation load for
% S6PPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% taug_reg [6x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PPRRRR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_gravloadJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:56:21
% EndTime: 2019-05-04 20:56:22
% DurationCPUTime: 0.31s
% Computational Cost: add. (512->79), mult. (1337->146), div. (0->0), fcn. (1763->16), ass. (0->53)
t51 = sin(pkin(13));
t52 = sin(pkin(12));
t42 = t52 * t51;
t55 = cos(pkin(13));
t56 = cos(pkin(12));
t49 = t56 * t55;
t58 = cos(pkin(6));
t35 = -t58 * t49 + t42;
t53 = sin(pkin(7));
t54 = sin(pkin(6));
t46 = t54 * t53;
t57 = cos(pkin(7));
t66 = t35 * t57 + t56 * t46;
t43 = t52 * t55;
t47 = t56 * t51;
t36 = t58 * t43 + t47;
t45 = t54 * t52;
t65 = t36 * t57 - t53 * t45;
t64 = t55 * t57 * t54 + t58 * t53;
t63 = cos(qJ(3));
t27 = qJ(5) + qJ(6);
t25 = sin(t27);
t32 = cos(qJ(4));
t62 = t25 * t32;
t26 = cos(t27);
t61 = t26 * t32;
t28 = sin(qJ(5));
t60 = t28 * t32;
t31 = cos(qJ(5));
t59 = t31 * t32;
t48 = t56 * t54;
t44 = t54 * t51;
t21 = -t58 * t42 + t49;
t30 = sin(qJ(3));
t10 = t21 * t63 - t65 * t30;
t14 = t64 * t30 + t63 * t44;
t15 = t35 * t53 - t57 * t48;
t16 = t36 * t53 + t57 * t45;
t19 = -t55 * t46 + t58 * t57;
t29 = sin(qJ(4));
t20 = t58 * t47 + t43;
t8 = t20 * t63 - t66 * t30;
t41 = g(1) * (-t10 * t29 + t16 * t32) + g(2) * (t15 * t32 - t8 * t29) + g(3) * (-t14 * t29 + t19 * t32);
t13 = t30 * t44 - t64 * t63;
t7 = t20 * t30 + t66 * t63;
t9 = t21 * t30 + t65 * t63;
t40 = g(1) * t9 + g(2) * t7 + g(3) * t13;
t12 = t14 * t32 + t19 * t29;
t6 = t10 * t32 + t16 * t29;
t4 = t15 * t29 + t8 * t32;
t2 = -g(1) * (-t9 * t25 - t6 * t26) - g(2) * (-t7 * t25 - t4 * t26) - g(3) * (-t12 * t26 - t13 * t25);
t1 = -g(1) * (-t6 * t25 + t9 * t26) - g(2) * (-t4 * t25 + t7 * t26) - g(3) * (-t12 * t25 + t13 * t26);
t3 = [-g(3), -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, -g(1) * t45 + g(2) * t48 - g(3) * t58, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, t40, g(1) * t10 + g(2) * t8 + g(3) * t14, 0, 0, 0, 0, 0, t40 * t32, -t40 * t29, 0, 0, 0, 0, 0, -g(1) * (t10 * t28 - t9 * t59) - g(2) * (t8 * t28 - t7 * t59) - g(3) * (-t13 * t59 + t14 * t28) -g(1) * (t10 * t31 + t9 * t60) - g(2) * (t8 * t31 + t7 * t60) - g(3) * (t13 * t60 + t14 * t31) 0, 0, 0, 0, 0, -g(1) * (t10 * t25 - t9 * t61) - g(2) * (t8 * t25 - t7 * t61) - g(3) * (-t13 * t61 + t14 * t25) -g(1) * (t10 * t26 + t9 * t62) - g(2) * (t8 * t26 + t7 * t62) - g(3) * (t13 * t62 + t14 * t26); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, g(1) * t6 + g(2) * t4 + g(3) * t12, 0, 0, 0, 0, 0, -t41 * t31, t41 * t28, 0, 0, 0, 0, 0, -t41 * t26, t41 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t6 * t28 + t9 * t31) - g(2) * (-t4 * t28 + t7 * t31) - g(3) * (-t12 * t28 + t13 * t31) -g(1) * (-t9 * t28 - t6 * t31) - g(2) * (-t7 * t28 - t4 * t31) - g(3) * (-t12 * t31 - t13 * t28) 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t3;
