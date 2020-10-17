% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% 
% Output:
% MM_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRPP3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:27:23
% EndTime: 2019-05-05 21:27:25
% DurationCPUTime: 0.58s
% Computational Cost: add. (461->102), mult. (810->173), div. (0->0), fcn. (762->6), ass. (0->66)
t46 = cos(qJ(4));
t32 = qJ(5) * t46;
t44 = sin(qJ(4));
t76 = -pkin(4) * t44 + t32;
t42 = pkin(4) + qJ(6);
t61 = t44 * qJ(5);
t75 = -t42 * t46 - t61;
t74 = -0.2e1 * t44;
t45 = sin(qJ(3));
t73 = 0.2e1 * t45;
t72 = pkin(3) * t46;
t47 = cos(qJ(3));
t70 = pkin(8) * t47;
t41 = cos(pkin(9));
t24 = -t41 * pkin(1) - pkin(2);
t13 = -t47 * pkin(3) - t45 * pkin(8) + t24;
t40 = sin(pkin(9));
t23 = t40 * pkin(1) + pkin(7);
t65 = t47 * t23;
t7 = t44 * t13 + t46 * t65;
t69 = t23 * t44;
t36 = t44 ^ 2;
t68 = t36 * t45;
t25 = t44 * t45;
t66 = t44 * t46;
t29 = t46 * t45;
t30 = t46 * t47;
t64 = -t46 * t13 + t44 * t65;
t21 = t45 * t32;
t63 = -pkin(4) * t25 + t21;
t38 = t46 ^ 2;
t62 = t36 + t38;
t60 = t47 * qJ(5);
t59 = t47 * t73;
t35 = t47 * pkin(4);
t4 = t35 + t64;
t18 = t45 * t23;
t8 = t18 - t63;
t58 = t62 * pkin(8);
t3 = t60 - t7;
t57 = -t3 * t46 + t4 * t44;
t55 = -t46 * pkin(4) - t61;
t17 = -pkin(3) + t55;
t56 = -t17 * t45 - t70;
t33 = t44 * pkin(8);
t19 = t44 * pkin(5) + t33;
t34 = t46 * pkin(8);
t20 = t46 * pkin(5) + t34;
t54 = t19 * t44 + t20 * t46;
t53 = -pkin(5) * t25 + t7;
t52 = -pkin(5) * t29 - t4;
t49 = qJ(5) ^ 2;
t48 = 0.2e1 * qJ(5);
t39 = t47 ^ 2;
t37 = t45 ^ 2;
t31 = -0.2e1 * t60;
t28 = t38 * t45;
t27 = t38 * t37;
t26 = t44 * t47;
t14 = t28 + t68;
t12 = -pkin(3) + t75;
t11 = t36 * t37 + t27 + t39;
t5 = qJ(6) * t25 + t8;
t2 = t53 - t60;
t1 = t47 * qJ(6) - t52;
t6 = [1, 0, 0 (t40 ^ 2 + t41 ^ 2) * pkin(1) ^ 2, t37, t59, 0, 0, 0, -0.2e1 * t24 * t47, t24 * t73, t27, -0.2e1 * t37 * t66, -0.2e1 * t45 * t30, t44 * t59, t39, 0.2e1 * t37 * t69 + 0.2e1 * t47 * t64, 0.2e1 * t37 * t23 * t46 + 0.2e1 * t7 * t47 (t3 * t44 + t4 * t46) * t73, -0.2e1 * t8 * t25 - 0.2e1 * t4 * t47, -0.2e1 * t8 * t29 + 0.2e1 * t3 * t47, t3 ^ 2 + t4 ^ 2 + t8 ^ 2 (t1 * t46 - t2 * t44) * t73, -0.2e1 * t2 * t47 - 0.2e1 * t5 * t29, 0.2e1 * t1 * t47 + 0.2e1 * t5 * t25, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57 * t45 - t8 * t47, 0, 0, 0, -t5 * t47 + (t1 * t44 + t2 * t46) * t45; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, 0, 0, t11; 0, 0, 0, 0, 0, 0, t45, t47, 0, -t18, -t65, t44 * t29, t28 - t68, -t26, -t30, 0, -t23 * t29 + (-pkin(3) * t45 + t70) * t44, pkin(8) * t30 + (t69 - t72) * t45, t57, t56 * t44 + t8 * t46, -t8 * t44 + t56 * t46, t57 * pkin(8) + t8 * t17 (t19 * t45 + t2) * t46 + (-t20 * t45 + t1) * t44, -t12 * t29 - t20 * t47 - t5 * t44, t12 * t25 + t19 * t47 - t5 * t46, t1 * t19 + t5 * t12 + t2 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, -t45, 0, 0, 0, 0, 0, t30, -t26, t14, -t30, t26, -t47 * t17 + t45 * t58, t14, t26, t30, -t47 * t12 + t54 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t36, 0.2e1 * t66, 0, 0, 0, 0.2e1 * t72, pkin(3) * t74, 0.2e1 * t58, 0.2e1 * t17 * t46, t17 * t74, t62 * pkin(8) ^ 2 + t17 ^ 2, 0.2e1 * t54, t12 * t74, -0.2e1 * t12 * t46, t12 ^ 2 + t19 ^ 2 + t20 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t25, -t47, -t64, -t7, t55 * t45, 0.2e1 * t35 + t64, t31 + t7, -t4 * pkin(4) - t3 * qJ(5), t75 * t45, t31 + t53 (-qJ(6) - t42) * t47 + t52, t2 * qJ(5) - t1 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, -t29, 0, t25, t29, t63, 0, t29, -t25, -t42 * t25 + t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t46, 0, -t33, -t34, t76, t33, t34, t76 * pkin(8), -t42 * t44 + t32, t20, -t19, t20 * qJ(5) - t19 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(4), t48, pkin(4) ^ 2 + t49, 0, t48, 0.2e1 * t42, t42 ^ 2 + t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t47, 0, t4, t29, 0, t47, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, 0, 0, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, 0, t33, t44, 0, 0, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(4), 0, 0, -1, -t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, -t47, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0, 0, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t6;
