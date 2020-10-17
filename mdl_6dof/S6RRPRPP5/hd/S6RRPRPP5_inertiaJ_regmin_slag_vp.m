% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRPP5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:44:48
% EndTime: 2019-05-06 12:44:51
% DurationCPUTime: 0.66s
% Computational Cost: add. (489->104), mult. (778->170), div. (0->0), fcn. (722->4), ass. (0->67)
t44 = sin(qJ(4));
t48 = pkin(4) + pkin(5);
t46 = cos(qJ(4));
t61 = t46 * qJ(5);
t77 = -t44 * t48 + t61;
t63 = t44 * qJ(5);
t20 = t48 * t46 + t63;
t45 = sin(qJ(2));
t76 = -0.2e1 * t45;
t47 = cos(qJ(2));
t75 = 0.2e1 * t47;
t74 = 0.2e1 * t48;
t73 = 2 * qJ(3);
t50 = -pkin(2) - pkin(8);
t72 = t45 * pkin(4);
t71 = t44 * t47;
t70 = t45 * t47;
t69 = t45 * t50;
t68 = t46 * t44;
t67 = t46 * t47;
t66 = t46 * t50;
t57 = -t45 * qJ(3) - pkin(1);
t15 = t50 * t47 + t57;
t36 = t45 * pkin(7);
t26 = t45 * pkin(3) + t36;
t65 = t44 * t15 - t46 * t26;
t8 = t46 * t15 + t44 * t26;
t37 = t47 * pkin(7);
t27 = t47 * pkin(3) + t37;
t39 = t44 ^ 2;
t41 = t46 ^ 2;
t29 = t39 + t41;
t40 = t45 ^ 2;
t42 = t47 ^ 2;
t64 = t40 + t42;
t62 = t44 * qJ(6);
t60 = t47 * qJ(3);
t59 = -0.2e1 * t70;
t35 = t45 * qJ(5);
t5 = t35 + t8;
t58 = 0.2e1 * t35 + t8;
t6 = t65 - t72;
t1 = t5 * t44 - t6 * t46;
t56 = -t45 * pkin(2) + t60;
t25 = pkin(4) * t46 + t63;
t55 = t44 * pkin(4) - t61;
t54 = -t47 * t62 - t65;
t52 = qJ(5) ^ 2;
t51 = 0.2e1 * qJ(5);
t33 = t44 * t50;
t32 = t46 * t45;
t31 = t44 * t45;
t30 = qJ(6) * t67;
t28 = t45 * t66;
t24 = -t47 * pkin(2) + t57;
t23 = qJ(3) + t55;
t21 = (-qJ(6) - t50) * t46;
t19 = t33 + t62;
t18 = t29 * t50;
t14 = -qJ(3) + t77;
t11 = t25 * t47 + t27;
t10 = t19 * t44 - t21 * t46;
t9 = -t20 * t47 - t27;
t4 = t30 + t5;
t3 = t4 * t44;
t2 = -t48 * t45 - t54;
t7 = [1, 0, 0, t40, 0.2e1 * t70, 0, 0, 0, pkin(1) * t75, pkin(1) * t76, 0.2e1 * t64 * pkin(7), t24 * t75, t24 * t76, t64 * pkin(7) ^ 2 + t24 ^ 2, t39 * t42, 0.2e1 * t42 * t68, t44 * t59, t46 * t59, t40, 0.2e1 * t27 * t67 - 0.2e1 * t45 * t65, -0.2e1 * t27 * t71 - 0.2e1 * t8 * t45, 0.2e1 * t11 * t67 - 0.2e1 * t6 * t45 (-t44 * t6 - t46 * t5) * t75, 0.2e1 * t11 * t71 + 0.2e1 * t5 * t45, t11 ^ 2 + t5 ^ 2 + t6 ^ 2, -0.2e1 * t2 * t45 - 0.2e1 * t9 * t67, 0.2e1 * t4 * t45 - 0.2e1 * t9 * t71 (t2 * t44 + t4 * t46) * t75, t2 ^ 2 + t4 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, t45, t47, 0, -t36, -t37, t56, t36, t37, t56 * pkin(7), -t44 * t67 (t39 - t41) * t47, t32, -t31, 0, t27 * t44 + t46 * t60 + t28, t27 * t46 + (-t60 - t69) * t44, t11 * t44 + t23 * t67 + t28, -t1, -t11 * t46 + (t23 * t47 + t69) * t44, t1 * t50 + t11 * t23, -t14 * t67 - t21 * t45 - t9 * t44, -t14 * t71 + t19 * t45 + t9 * t46, t21 * t71 + t3 + (t19 * t47 - t2) * t46, t9 * t14 + t4 * t19 + t2 * t21; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(2), t73, pkin(2) ^ 2 + (qJ(3) ^ 2) t41, -0.2e1 * t68, 0, 0, 0, t44 * t73, t46 * t73, 0.2e1 * t23 * t44, -0.2e1 * t18, -0.2e1 * t23 * t46, t29 * t50 ^ 2 + t23 ^ 2, -0.2e1 * t14 * t44, 0.2e1 * t14 * t46, 0.2e1 * t10, t14 ^ 2 + t19 ^ 2 + t21 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, 0, t36, 0, 0, 0, 0, 0, t32, -t31, t32, 0, t31, t1, t32, t31, 0, -t2 * t46 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, -t29, 0, t18, 0, 0, t29, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, 0, 0, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, -t67, t45, -t65, -t8, -t65 + 0.2e1 * t72, t55 * t47, t58, -t6 * pkin(4) + t5 * qJ(5), t45 * t74 + t54, t30 + t58, t77 * t47, t4 * qJ(5) - t2 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t44, 0, t66, -t33, t66, -t25, t33, t25 * t50, -t21, t19, t20, t19 * qJ(5) - t21 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t44, t46, 0, t44, t25, t46, t44, 0, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, t51, pkin(4) ^ 2 + t52, t74, t51, 0, t48 ^ 2 + t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, -t71, 0, t6, -t45, 0, t71, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0, -t66, 0, 0, -t46, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, 0, 0, 0, -t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4), -1, 0, 0, -t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, -t71, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, t46, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t7;
