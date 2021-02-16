% Calculate minimal parameter regressor of gravitation load for
% S5RRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% taug_reg [5x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:38
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPP6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:37:07
% EndTime: 2021-01-15 22:37:09
% DurationCPUTime: 0.27s
% Computational Cost: add. (233->69), mult. (372->103), div. (0->0), fcn. (387->10), ass. (0->57)
t41 = sin(pkin(8));
t42 = cos(pkin(8));
t27 = pkin(4) * t42 + qJ(5) * t41 + pkin(3);
t28 = -t41 * pkin(4) + qJ(5) * t42;
t44 = sin(qJ(3));
t47 = cos(qJ(3));
t9 = -t27 * t44 + t28 * t47;
t71 = -pkin(6) + t9;
t70 = -t27 * t47 - t28 * t44;
t49 = cos(qJ(1));
t56 = t49 * t47;
t46 = sin(qJ(1));
t48 = cos(qJ(2));
t61 = t48 * t46;
t18 = t44 * t61 + t56;
t57 = t49 * t44;
t20 = t46 * t47 - t48 * t57;
t45 = sin(qJ(2));
t66 = g(3) * t45;
t69 = -g(1) * t20 + g(2) * t18 + t44 * t66;
t52 = g(1) * t49 + g(2) * t46;
t14 = -g(3) * t48 + t52 * t45;
t43 = qJ(4) + pkin(7);
t32 = t43 * t45;
t64 = pkin(1) + t32;
t60 = t48 * t49;
t40 = qJ(3) + pkin(8);
t36 = sin(t40);
t59 = t49 * t36;
t37 = cos(t40);
t58 = t49 * t37;
t10 = t36 * t61 + t58;
t12 = -t46 * t37 + t48 * t59;
t53 = g(1) * t10 - g(2) * t12;
t51 = g(1) * t46 - g(2) * t49;
t11 = t37 * t61 - t59;
t13 = t46 * t36 + t48 * t58;
t50 = g(1) * t13 + g(2) * t11 + t37 * t66;
t35 = t47 * pkin(3) + pkin(2);
t34 = t44 * pkin(3) + pkin(6);
t33 = t43 * t48;
t29 = t35 * t48;
t26 = t41 * t47 + t42 * t44;
t25 = t41 * t44 - t42 * t47;
t22 = t51 * t45;
t21 = t46 * t44 + t48 * t56;
t19 = -t47 * t61 + t57;
t16 = t29 + t64;
t15 = t52 * t48 + t66;
t8 = pkin(2) - t70;
t7 = t8 * t48;
t6 = t14 * t37;
t5 = t14 * t36;
t3 = t7 + t64;
t2 = g(1) * t11 - g(2) * t13;
t1 = g(1) * t12 + g(2) * t10 + t36 * t66;
t4 = [0, t51, t52, 0, 0, 0, 0, 0, t51 * t48, -t22, 0, 0, 0, 0, 0, -g(1) * t19 - g(2) * t21, -g(1) * t18 - g(2) * t20, t2, -t53, t22, -g(1) * (-t16 * t46 + t34 * t49) - g(2) * (t16 * t49 + t34 * t46), t2, t22, t53, -g(1) * (-t3 * t46 - t71 * t49) - g(2) * (t3 * t49 - t71 * t46); 0, 0, 0, 0, 0, 0, 0, 0, t14, t15, 0, 0, 0, 0, 0, t14 * t47, -t14 * t44, t6, -t5, -t15, -g(3) * (t29 + t32) - t52 * (-t45 * t35 + t33), t6, -t15, t5, -g(3) * (t7 + t32) - t52 * (-t8 * t45 + t33); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, g(1) * t21 - g(2) * t19 + t47 * t66, t1, t50, 0, t69 * pkin(3), t1, 0, -t50, -g(1) * (-t46 * t70 + t9 * t60) - g(2) * (t70 * t49 + t9 * t61) - t9 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t25 * t46 + t26 * t60) - g(2) * (-t25 * t49 + t26 * t61) - t26 * t66;];
taug_reg = t4;
