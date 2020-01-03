% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x23]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRP7_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP7_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:01:31
% EndTime: 2019-12-31 20:01:34
% DurationCPUTime: 0.69s
% Computational Cost: add. (1157->120), mult. (2604->243), div. (0->0), fcn. (2310->6), ass. (0->82)
t45 = sin(pkin(8));
t46 = cos(pkin(8));
t48 = sin(qJ(2));
t50 = cos(qJ(2));
t31 = t45 * t48 - t46 * t50;
t32 = t45 * t50 + t46 * t48;
t70 = -pkin(2) * t50 - pkin(1);
t20 = pkin(3) * t31 - pkin(7) * t32 + t70;
t85 = -qJ(3) - pkin(6);
t34 = t85 * t48;
t35 = t85 * t50;
t23 = t34 * t45 - t35 * t46;
t47 = sin(qJ(4));
t49 = cos(qJ(4));
t99 = t47 * t20 + t49 * t23;
t67 = qJD(2) * t85;
t25 = t50 * qJD(3) + t48 * t67;
t55 = -t48 * qJD(3) + t50 * t67;
t13 = t46 * t25 + t45 * t55;
t27 = t32 * qJD(2);
t76 = t50 * qJD(2);
t78 = t48 * qJD(2);
t28 = -t45 * t78 + t46 * t76;
t41 = pkin(2) * t78;
t14 = pkin(3) * t27 - pkin(7) * t28 + t41;
t4 = -qJD(4) * t99 - t47 * t13 + t14 * t49;
t62 = pkin(4) * t49 + qJ(5) * t47;
t98 = qJD(4) * t62 - t49 * qJD(5);
t42 = qJD(4) * t49;
t81 = qJD(4) * t47;
t3 = -t13 * t49 - t14 * t47 - t20 * t42 + t23 * t81;
t80 = t31 * qJD(5);
t82 = t27 * qJ(5);
t1 = -t3 + t80 + t82;
t94 = t27 * pkin(4);
t2 = -t4 - t94;
t6 = qJ(5) * t31 + t99;
t60 = t20 * t49 - t23 * t47;
t7 = -pkin(4) * t31 - t60;
t97 = t1 * t47 - t2 * t49 + (t47 * t7 + t49 * t6) * qJD(4);
t96 = 0.2e1 * qJD(4);
t95 = 0.2e1 * qJD(5);
t12 = t45 * t25 - t46 * t55;
t93 = t12 * t47;
t38 = pkin(2) * t45 + pkin(7);
t92 = t27 * t38;
t91 = t31 * t38;
t90 = t32 * t49;
t89 = t47 * t27;
t88 = t47 * t28;
t87 = t49 * t27;
t86 = t49 * t28;
t43 = t47 ^ 2;
t44 = t49 ^ 2;
t83 = t43 - t44;
t79 = t47 * qJD(5);
t75 = -0.2e1 * pkin(1) * qJD(2);
t39 = -pkin(2) * t46 - pkin(3);
t74 = t39 * t96;
t73 = t38 * t81;
t72 = t38 * t42;
t71 = t47 * t42;
t69 = -0.4e1 * t47 * t90;
t22 = -t46 * t34 - t45 * t35;
t66 = t83 * qJD(4);
t63 = -t47 * t6 + t49 * t7;
t61 = pkin(4) * t47 - qJ(5) * t49;
t58 = t28 * t39 - t92;
t57 = -t32 * t39 + t91;
t19 = t31 * t42 + t89;
t17 = -t31 * t81 + t87;
t56 = t32 * t42 + t88;
t18 = -t32 * t81 + t86;
t29 = t39 - t62;
t5 = t61 * t28 + t32 * t98 + t12;
t54 = -t5 + (t29 * t32 - t91) * qJD(4);
t26 = -pkin(4) * t81 + qJ(5) * t42 + t79;
t8 = t32 * t61 + t22;
t53 = qJD(4) * t8 - t26 * t32 + t28 * t29 - t92;
t51 = qJD(4) * t63 + t1 * t49 + t2 * t47;
t30 = t32 ^ 2;
t9 = [0, 0, 0, 0.2e1 * t48 * t76, 0.2e1 * (-t48 ^ 2 + t50 ^ 2) * qJD(2), 0, 0, 0, t48 * t75, t50 * t75, 0.2e1 * t12 * t32 - 0.2e1 * t13 * t31 + 0.2e1 * t22 * t28 - 0.2e1 * t23 * t27, 0.2e1 * t22 * t12 + 0.2e1 * t23 * t13 + 0.2e1 * t41 * t70, 0.2e1 * t28 * t32 * t44 - 0.2e1 * t30 * t71, t30 * t83 * t96 + t28 * t69, 0.2e1 * t18 * t31 + 0.2e1 * t32 * t87, -0.2e1 * t31 * t56 - 0.2e1 * t32 * t89, 0.2e1 * t31 * t27, 0.2e1 * t22 * t56 + 0.2e1 * t27 * t60 + 0.2e1 * t4 * t31 + 0.2e1 * t32 * t93, 0.2e1 * t12 * t90 + 0.2e1 * t18 * t22 - 0.2e1 * t27 * t99 + 0.2e1 * t3 * t31, 0.2e1 * t8 * t88 - 0.2e1 * t2 * t31 - 0.2e1 * t7 * t27 + 0.2e1 * (t42 * t8 + t5 * t47) * t32, 0.2e1 * t63 * t28 - 0.2e1 * t32 * t97, -0.2e1 * t8 * t86 + 0.2e1 * t1 * t31 + 0.2e1 * t6 * t27 + 0.2e1 * (-t5 * t49 + t8 * t81) * t32, 0.2e1 * t1 * t6 + 0.2e1 * t2 * t7 + 0.2e1 * t5 * t8; 0, 0, 0, 0, 0, t76, -t78, 0, -pkin(6) * t76, pkin(6) * t78, (-t27 * t45 - t28 * t46) * pkin(2), (-t12 * t46 + t13 * t45) * pkin(2), -t32 * t66 + t47 * t86, qJD(4) * t69 - t28 * t83, t19, t17, 0, -t12 * t49 + t58 * t47 + (t22 * t47 - t49 * t57) * qJD(4), t93 + t58 * t49 + (t22 * t49 + t47 * t57) * qJD(4), t47 * t53 + t49 * t54, t51, t47 * t54 - t49 * t53, -t8 * t26 + t5 * t29 + t38 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t71, -0.2e1 * t66, 0, 0, 0, t47 * t74, t49 * t74, 0.2e1 * t26 * t49 + 0.2e1 * t29 * t81, 0, 0.2e1 * t26 * t47 - 0.2e1 * t29 * t42, -0.2e1 * t29 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, 0, 0, 0, 0, t17, -t19, t17, (-t43 - t44) * t28, t19, t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t56, t27, t4, t3, t4 + 0.2e1 * t94, -t62 * t28 + (qJD(4) * t61 - t79) * t32, -t3 + 0.2e1 * t80 + 0.2e1 * t82, -pkin(4) * t2 + qJ(5) * t1 + qJD(5) * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, -t81, 0, -t72, t73, -t72, -t98, -t73, -t98 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, -t42, -t81, 0, t42, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, qJ(5) * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, t18, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;
