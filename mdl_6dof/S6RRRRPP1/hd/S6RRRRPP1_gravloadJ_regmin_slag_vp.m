% Calculate minimal parameter regressor of gravitation load for
% S6RRRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
% 
% Output:
% taug_reg [6x30]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP1_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 17:57:05
% EndTime: 2019-05-07 17:57:07
% DurationCPUTime: 0.36s
% Computational Cost: add. (412->80), mult. (450->109), div. (0->0), fcn. (454->10), ass. (0->58)
t42 = cos(qJ(2));
t41 = cos(qJ(4));
t28 = t41 * pkin(4) + pkin(3);
t36 = qJ(2) + qJ(3);
t33 = cos(t36);
t19 = t33 * t28;
t32 = sin(t36);
t37 = -qJ(5) - pkin(9);
t56 = -t32 * t37 + t19;
t55 = t42 * pkin(2) + t56;
t35 = qJ(4) + pkin(10);
t30 = sin(t35);
t31 = cos(t35);
t77 = pkin(5) * t31 + qJ(6) * t30;
t76 = -pkin(1) - t55;
t75 = t77 * t33;
t40 = sin(qJ(1));
t43 = cos(qJ(1));
t54 = g(1) * t43 + g(2) * t40;
t5 = -g(3) * t33 + t54 * t32;
t39 = sin(qJ(2));
t74 = pkin(2) * t39;
t70 = g(3) * t32;
t67 = t33 * t37;
t66 = t40 * t30;
t65 = t40 * t31;
t38 = sin(qJ(4));
t64 = t40 * t38;
t63 = t40 * t41;
t62 = t43 * t30;
t61 = t43 * t31;
t60 = t43 * t38;
t59 = t43 * t41;
t57 = t33 * t60;
t53 = g(1) * t40 - g(2) * t43;
t52 = t28 * t32 + t67;
t12 = t33 * t64 + t59;
t51 = t28 + t77;
t50 = t54 * t33;
t7 = t33 * t66 + t61;
t9 = t33 * t62 - t65;
t48 = g(1) * t9 + g(2) * t7 + t30 * t70;
t44 = -pkin(8) - pkin(7);
t47 = pkin(4) * t64 - t40 * t44 - t76 * t43;
t6 = t50 + t70;
t46 = pkin(4) * t60 + t76 * t40 - t43 * t44;
t25 = pkin(4) * t63;
t15 = t33 * t59 + t64;
t14 = -t57 + t63;
t13 = -t33 * t63 + t60;
t11 = t53 * t32;
t10 = t33 * t61 + t66;
t8 = t33 * t65 - t62;
t4 = t5 * t41;
t3 = t5 * t38;
t2 = t5 * t31;
t1 = t5 * t30;
t16 = [0, t53, t54, 0, 0, 0, 0, 0, t53 * t42, -t53 * t39, 0, 0, 0, 0, 0, t53 * t33, -t11, 0, 0, 0, 0, 0, -g(1) * t13 - g(2) * t15, -g(1) * t12 - g(2) * t14, t11, -g(1) * t46 - g(2) * t47, g(1) * t8 - g(2) * t10, t11, g(1) * t7 - g(2) * t9, -g(1) * (-t8 * pkin(5) - t7 * qJ(6) + t46) - g(2) * (t10 * pkin(5) + t9 * qJ(6) + t47); 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t42 + t54 * t39, g(3) * t39 + t54 * t42, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(3) * t55 + t54 * (t52 + t74) t2, -t6, t1, -g(3) * (t55 + t75) + t54 * (t51 * t32 + t67 + t74); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(3) * t56 + t54 * t52, t2, -t6, t1, -g(3) * (t19 + t75) + t37 * t50 + (g(3) * t37 + t54 * t51) * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t14 + g(2) * t12 + t38 * t70, g(1) * t15 - g(2) * t13 + t41 * t70, 0, -g(1) * t25 + (g(2) * t59 + t6 * t38) * pkin(4), t48, 0, -g(1) * t10 - g(2) * t8 - t31 * t70, -g(1) * (-pkin(4) * t57 - t9 * pkin(5) + t10 * qJ(6) + t25) - g(2) * (-t12 * pkin(4) - t7 * pkin(5) + t8 * qJ(6)) - (-pkin(4) * t38 - pkin(5) * t30 + qJ(6) * t31) * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48;];
taug_reg  = t16;
