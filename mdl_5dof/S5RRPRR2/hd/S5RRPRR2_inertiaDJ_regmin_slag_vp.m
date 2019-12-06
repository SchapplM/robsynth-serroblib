% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x26]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRR2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR2_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:28:26
% EndTime: 2019-12-05 18:28:28
% DurationCPUTime: 0.48s
% Computational Cost: add. (1345->86), mult. (3010->167), div. (0->0), fcn. (2983->8), ass. (0->70)
t61 = sin(qJ(4));
t64 = cos(qJ(4));
t59 = cos(pkin(9));
t74 = t59 * pkin(2) + pkin(3);
t58 = sin(pkin(9));
t86 = pkin(2) * t58;
t87 = t61 * t86 - t64 * t74;
t62 = sin(qJ(2));
t65 = cos(qJ(2));
t83 = -qJ(3) - pkin(6);
t71 = qJD(2) * t83;
t43 = t65 * qJD(3) + t62 * t71;
t44 = -t62 * qJD(3) + t65 * t71;
t23 = -t58 * t43 + t59 * t44;
t80 = t65 * qJD(2);
t81 = t62 * qJD(2);
t47 = -t58 * t81 + t59 * t80;
t20 = -t47 * pkin(7) + t23;
t24 = t59 * t43 + t58 * t44;
t49 = t58 * t65 + t59 * t62;
t46 = t49 * qJD(2);
t21 = -t46 * pkin(7) + t24;
t54 = t83 * t62;
t55 = t83 * t65;
t29 = t59 * t54 + t58 * t55;
t25 = -t49 * pkin(7) + t29;
t30 = t58 * t54 - t59 * t55;
t48 = -t58 * t62 + t59 * t65;
t26 = t48 * pkin(7) + t30;
t68 = t64 * t25 - t61 * t26;
t7 = -t68 * qJD(4) - t61 * t20 - t64 * t21;
t41 = pkin(4) - t87;
t85 = -pkin(4) - t41;
t42 = t61 * t74 + t64 * t86;
t63 = cos(qJ(5));
t84 = t42 * t63;
t60 = sin(qJ(5));
t82 = qJD(5) * t60;
t78 = -0.2e1 * pkin(1) * qJD(2);
t57 = pkin(2) * t81;
t77 = pkin(4) * t82;
t76 = qJD(5) * t63 * pkin(4);
t75 = -t65 * pkin(2) - pkin(1);
t31 = t46 * pkin(3) + t57;
t38 = t87 * qJD(4);
t39 = t42 * qJD(4);
t73 = t60 * t38 - t63 * t39;
t72 = t60 * t39 + t42 * t82;
t67 = -t61 * t25 - t64 * t26;
t28 = t61 * t48 + t64 * t49;
t66 = t64 * t48 - t61 * t49;
t17 = t60 * t28 - t63 * t66;
t18 = t63 * t28 + t60 * t66;
t34 = -t48 * pkin(3) + t75;
t8 = t67 * qJD(4) + t64 * t20 - t61 * t21;
t22 = -pkin(4) * t66 + t34;
t16 = t28 * qJD(4) + t64 * t46 + t61 * t47;
t15 = t66 * qJD(4) - t61 * t46 + t64 * t47;
t13 = (-t41 * t60 - t84) * qJD(5) + t73;
t12 = (-qJD(5) * t41 + t38) * t63 + t72;
t11 = t16 * pkin(4) + t31;
t10 = pkin(8) * t66 - t67;
t9 = -t28 * pkin(8) + t68;
t6 = t18 * qJD(5) + t60 * t15 + t63 * t16;
t5 = -t17 * qJD(5) + t63 * t15 - t60 * t16;
t4 = -t15 * pkin(8) + t8;
t3 = -t16 * pkin(8) - t7;
t2 = -t60 * t3 + t63 * t4 + (-t10 * t63 - t60 * t9) * qJD(5);
t1 = -t63 * t3 - t60 * t4 + (t10 * t60 - t63 * t9) * qJD(5);
t14 = [0, 0, 0, 0.2e1 * t62 * t80, 0.2e1 * (-t62 ^ 2 + t65 ^ 2) * qJD(2), 0, 0, 0, t62 * t78, t65 * t78, -0.2e1 * t23 * t49 + 0.2e1 * t24 * t48 - 0.2e1 * t29 * t47 - 0.2e1 * t30 * t46, 0.2e1 * t29 * t23 + 0.2e1 * t30 * t24 + 0.2e1 * t75 * t57, 0.2e1 * t28 * t15, 0.2e1 * t15 * t66 - 0.2e1 * t28 * t16, 0, 0, 0, 0.2e1 * t34 * t16 - 0.2e1 * t31 * t66, 0.2e1 * t34 * t15 + 0.2e1 * t31 * t28, 0.2e1 * t18 * t5, -0.2e1 * t5 * t17 - 0.2e1 * t18 * t6, 0, 0, 0, 0.2e1 * t11 * t17 + 0.2e1 * t22 * t6, 0.2e1 * t11 * t18 + 0.2e1 * t22 * t5; 0, 0, 0, 0, 0, t80, -t81, 0, -pkin(6) * t80, pkin(6) * t81, (-t46 * t58 - t47 * t59) * pkin(2), (t23 * t59 + t24 * t58) * pkin(2), 0, 0, t15, -t16, 0, t8, t7, 0, 0, t5, -t6, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t39, 0.2e1 * t38, 0, 0, 0, 0, 0, 0.2e1 * t13, 0.2e1 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, 0, 0, 0, 0, 0, t16, t15, 0, 0, 0, 0, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t16, 0, t8, t7, 0, 0, t5, -t6, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, t38, 0, 0, 0, 0, 0, (t85 * t60 - t84) * qJD(5) + t73, (t85 * qJD(5) + t38) * t63 + t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t77, -0.2e1 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t6, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77, -t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t14;
