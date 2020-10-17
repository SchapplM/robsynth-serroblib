% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRPRP3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:39:20
% EndTime: 2019-05-05 17:39:22
% DurationCPUTime: 0.60s
% Computational Cost: add. (712->106), mult. (1360->182), div. (0->0), fcn. (1458->8), ass. (0->63)
t40 = sin(pkin(10));
t42 = cos(pkin(10));
t44 = sin(qJ(5));
t70 = cos(qJ(5));
t78 = -t44 * t40 + t70 * t42;
t26 = t70 * t40 + t44 * t42;
t77 = -0.2e1 * t26;
t35 = -t42 * pkin(4) - pkin(3);
t76 = 0.2e1 * t35;
t45 = sin(qJ(3));
t75 = 0.2e1 * t45;
t46 = cos(qJ(3));
t74 = -0.2e1 * t46;
t73 = 0.2e1 * t46;
t43 = cos(pkin(9));
t34 = -t43 * pkin(1) - pkin(2);
t72 = t46 * pkin(3);
t23 = -t45 * qJ(4) + t34 - t72;
t41 = sin(pkin(9));
t33 = t41 * pkin(1) + pkin(7);
t58 = t46 * t42;
t14 = t40 * t23 + t33 * t58;
t64 = t40 * t45;
t10 = -pkin(8) * t64 + t14;
t21 = t42 * t23;
t63 = t42 * t45;
t65 = t33 * t40;
t8 = -pkin(8) * t63 + t21 + (-pkin(4) - t65) * t46;
t4 = t70 * t10 + t44 * t8;
t71 = t46 * pkin(5);
t57 = pkin(8) + qJ(4);
t28 = t57 * t42;
t53 = t57 * t40;
t15 = t44 * t28 + t70 * t53;
t69 = t15 * t46;
t16 = t70 * t28 - t44 * t53;
t68 = t16 * t46;
t18 = t26 * t45;
t67 = t18 * t26;
t66 = t26 * t46;
t29 = t45 * t33;
t61 = t46 * t78;
t60 = t46 * t33;
t59 = t46 * t40;
t22 = pkin(4) * t64 + t29;
t56 = t40 ^ 2 + t42 ^ 2;
t55 = t46 * qJ(6);
t3 = -t44 * t10 + t70 * t8;
t52 = t56 * qJ(4);
t51 = -pkin(3) * t45 + qJ(4) * t46;
t19 = t78 * t45;
t50 = -t18 * pkin(5) + t19 * qJ(6);
t13 = -t33 * t59 + t21;
t49 = -t13 * t40 + t14 * t42;
t39 = t46 ^ 2;
t38 = t45 ^ 2;
t17 = t19 ^ 2;
t12 = t19 * t78;
t11 = -pkin(5) * t78 - t26 * qJ(6) + t35;
t5 = -t50 + t22;
t2 = -t3 + t71;
t1 = -t55 + t4;
t6 = [1, 0, 0 (t41 ^ 2 + t43 ^ 2) * pkin(1) ^ 2, t38, t45 * t73, 0, 0, 0, t34 * t74, t34 * t75, -0.2e1 * t13 * t46 + 0.2e1 * t38 * t65, 0.2e1 * t38 * t33 * t42 + 0.2e1 * t14 * t46 (-t13 * t42 - t14 * t40) * t75, t38 * t33 ^ 2 + t13 ^ 2 + t14 ^ 2, t17, -0.2e1 * t18 * t19, t19 * t74, t18 * t73, t39, 0.2e1 * t22 * t18 - 0.2e1 * t3 * t46, 0.2e1 * t22 * t19 + 0.2e1 * t4 * t46, 0.2e1 * t5 * t18 + 0.2e1 * t2 * t46, -0.2e1 * t1 * t18 + 0.2e1 * t2 * t19, -0.2e1 * t1 * t46 - 0.2e1 * t5 * t19, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (t49 - t60) * t45, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 * t19 + t2 * t18 - t5 * t46; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56 * t38 + t39, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18 ^ 2 + t17 + t39; 0, 0, 0, 0, 0, 0, t45, t46, 0, -t29, -t60, -t42 * t29 + t51 * t40, t40 * t29 + t51 * t42, t49, -pkin(3) * t29 + t49 * qJ(4), t19 * t26, t12 - t67, -t66, -t61, 0, t35 * t18 - t22 * t78 + t69, t35 * t19 + t22 * t26 + t68, t11 * t18 - t5 * t78 + t69, t1 * t78 + t15 * t19 - t16 * t18 + t2 * t26, -t11 * t19 - t5 * t26 - t68, t1 * t16 + t5 * t11 + t2 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t45, t58, -t59, t56 * t45, t45 * t52 + t72, 0, 0, 0, 0, 0, t61, -t66, t61, t12 + t67, t66, -t46 * t11 + t18 * t15 + t19 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3) * t42, -0.2e1 * pkin(3) * t40, 0.2e1 * t52, t56 * qJ(4) ^ 2 + pkin(3) ^ 2, t26 ^ 2, -t78 * t77, 0, 0, 0, -t78 * t76, t26 * t76, -0.2e1 * t11 * t78, 0.2e1 * t15 * t26 + 0.2e1 * t16 * t78, t11 * t77, t11 ^ 2 + t15 ^ 2 + t16 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, t63, 0, t29, 0, 0, 0, 0, 0, t18, t19, t18, 0, -t19, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, t40, 0, -pkin(3), 0, 0, 0, 0, 0, -t78, t26, -t78, 0, -t26, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t18, -t46, t3, -t4, t3 - 0.2e1 * t71, -pkin(5) * t19 - t18 * qJ(6), -0.2e1 * t55 + t4, -t2 * pkin(5) + t1 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, -t19, -t18, 0, t19, t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, t78, 0, -t15, -t16, -t15, -pkin(5) * t26 + qJ(6) * t78, t16, -t15 * pkin(5) + t16 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, t19, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t6;
