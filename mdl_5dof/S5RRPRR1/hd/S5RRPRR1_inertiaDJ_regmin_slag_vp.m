% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:13
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR1_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_inertiaDJ_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:13:18
% EndTime: 2021-01-15 21:13:21
% DurationCPUTime: 0.52s
% Computational Cost: add. (516->87), mult. (1283->175), div. (0->0), fcn. (1085->6), ass. (0->78)
t45 = sin(qJ(2));
t40 = t45 ^ 2;
t48 = cos(qJ(2));
t42 = t48 ^ 2;
t87 = (t40 - t42) * qJD(2);
t46 = cos(qJ(5));
t41 = t46 ^ 2;
t43 = sin(qJ(5));
t77 = t43 ^ 2 - t41;
t61 = t77 * qJD(5);
t86 = qJD(2) + qJD(4);
t49 = pkin(1) + pkin(2);
t78 = pkin(3) + qJ(3);
t29 = t78 * t45;
t30 = t78 * t48;
t44 = sin(qJ(4));
t47 = cos(qJ(4));
t16 = -t44 * t29 + t47 * t30;
t62 = qJD(2) * t78;
t69 = t48 * qJD(3);
t51 = -t45 * t62 + t69;
t70 = t45 * qJD(3);
t52 = t48 * t62 + t70;
t7 = t16 * qJD(4) + t44 * t51 + t47 * t52;
t4 = t7 * t43;
t79 = t47 * t29;
t15 = t44 * t30 + t79;
t38 = qJD(5) * t46;
t85 = t15 * t38 + t4;
t26 = t44 * t45 - t47 * t48;
t13 = t86 * t26;
t27 = t44 * t48 + t47 * t45;
t84 = t27 * t13;
t83 = t27 * t46;
t14 = t86 * t27;
t82 = t43 * t14;
t81 = t46 * t13;
t80 = t46 * t14;
t71 = t45 * qJD(2);
t36 = pkin(1) * t71;
t28 = pkin(2) * t71 + t36;
t75 = qJD(4) * t44;
t74 = qJD(4) * t49;
t73 = qJD(5) * t43;
t72 = qJD(5) * t47;
t37 = t48 * qJD(2);
t68 = qJ(3) * qJD(2);
t31 = t49 * t48;
t67 = t44 * t74;
t66 = t47 * t74;
t65 = t43 * t38;
t64 = t45 * t37;
t63 = -0.4e1 * t43 * t83;
t60 = (t40 + t42) * qJD(3);
t17 = -t27 * pkin(4) - t31;
t59 = t46 * t16 + t43 * t17;
t58 = t43 * t16 - t46 * t17;
t33 = t44 * t49 + pkin(4);
t57 = t27 * t47 * t49 + t26 * t33;
t56 = -t43 * t13 + t27 * t38;
t55 = -t27 * t73 - t81;
t9 = t26 * t38 + t82;
t54 = t26 * t73 - t80;
t53 = -t48 * t68 - t70;
t50 = -t14 * t33 + (t13 * t47 + (-t26 * t47 + t27 * t44) * qJD(4)) * t49;
t32 = 0.2e1 * t65;
t25 = -0.2e1 * t61;
t24 = t27 ^ 2;
t19 = (-t43 * t72 - t46 * t75) * t49;
t18 = (t43 * t75 - t46 * t72) * t49;
t11 = t15 * t73;
t10 = t13 * pkin(4) + t28;
t6 = qJD(4) * t79 + t30 * t75 + t44 * t52 - t47 * t51;
t5 = -t27 * t61 - t43 * t81;
t3 = qJD(5) * t63 + t77 * t13;
t2 = -t59 * qJD(5) + t46 * t10 + t43 * t6;
t1 = t58 * qJD(5) - t43 * t10 + t46 * t6;
t8 = [0, 0, 0, 0.2e1 * t64, -0.2e1 * t87, 0, 0, 0, 0, 0, -0.4e1 * pkin(1) * t64, 0.2e1 * pkin(1) * t87, 0.2e1 * t60, -0.2e1 * pkin(1) ^ 2 * t64 + 0.2e1 * qJ(3) * t60, -0.2e1 * t84, 0.2e1 * t13 * t26 - 0.2e1 * t27 * t14, 0, 0, 0, -0.2e1 * t31 * t14 + 0.2e1 * t28 * t26, 0.2e1 * t31 * t13 + 0.2e1 * t28 * t27, -0.2e1 * t24 * t65 - 0.2e1 * t41 * t84, -t13 * t63 + 0.2e1 * t24 * t61, 0.2e1 * t55 * t26 + 0.2e1 * t27 * t80, -0.2e1 * t56 * t26 - 0.2e1 * t27 * t82, 0.2e1 * t26 * t14, -0.2e1 * t58 * t14 + 0.2e1 * t56 * t15 + 0.2e1 * t2 * t26 + 0.2e1 * t27 * t4, 0.2e1 * t1 * t26 - 0.2e1 * t59 * t14 + 0.2e1 * t55 * t15 + 0.2e1 * t7 * t83; 0, 0, 0, 0, 0, t37, -t71, 0, 0, 0, t53, t45 * t68 - t69, -pkin(1) * t37, t53 * pkin(1), 0, 0, -t13, -t14, 0, -t7, t6, t5, t3, t9, -t54, 0, t11 + (-t57 * qJD(5) - t7) * t46 + t50 * t43, t50 * t46 + t57 * t73 + t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t67, -0.2e1 * t66, t32, t25, 0, 0, 0, 0.2e1 * t19, 0.2e1 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, t37, 0, t36, 0, 0, 0, 0, 0, t14, -t13, 0, 0, 0, 0, 0, -t54, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t14, 0, -t7, t6, t5, t3, t9, -t54, 0, -t9 * pkin(4) - t7 * t46 + t11, t54 * pkin(4) + t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, -t66, t32, t25, 0, 0, 0, t19, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, t25, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, -t56, t14, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t73, 0, -t33 * t38 - t43 * t66, t33 * t73 - t46 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73, -t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t73, 0, -pkin(4) * t38, pkin(4) * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t8;
