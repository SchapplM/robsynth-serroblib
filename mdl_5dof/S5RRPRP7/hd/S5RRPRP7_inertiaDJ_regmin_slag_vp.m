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
% MMD_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:42
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
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
% StartTime: 2021-01-15 20:41:23
% EndTime: 2021-01-15 20:41:29
% DurationCPUTime: 0.74s
% Computational Cost: add. (1177->124), mult. (2662->251), div. (0->0), fcn. (2354->6), ass. (0->82)
t46 = sin(pkin(8));
t47 = cos(pkin(8));
t49 = sin(qJ(2));
t51 = cos(qJ(2));
t31 = t46 * t49 - t47 * t51;
t32 = t46 * t51 + t47 * t49;
t40 = -pkin(2) * t51 - pkin(1);
t20 = pkin(3) * t31 - pkin(7) * t32 + t40;
t85 = -qJ(3) - pkin(6);
t34 = t85 * t49;
t35 = t85 * t51;
t23 = t34 * t46 - t35 * t47;
t48 = sin(qJ(4));
t50 = cos(qJ(4));
t99 = t48 * t20 + t50 * t23;
t68 = qJD(2) * t85;
t25 = t51 * qJD(3) + t49 * t68;
t56 = -t49 * qJD(3) + t51 * t68;
t13 = t47 * t25 + t46 * t56;
t27 = t32 * qJD(2);
t76 = t51 * qJD(2);
t78 = t49 * qJD(2);
t28 = -t46 * t78 + t47 * t76;
t42 = pkin(2) * t78;
t14 = pkin(3) * t27 - pkin(7) * t28 + t42;
t4 = -qJD(4) * t99 - t48 * t13 + t14 * t50;
t63 = pkin(4) * t50 + qJ(5) * t48;
t98 = qJD(4) * t63 - t50 * qJD(5);
t43 = qJD(4) * t50;
t81 = qJD(4) * t48;
t3 = -t13 * t50 - t14 * t48 - t20 * t43 + t23 * t81;
t80 = t31 * qJD(5);
t82 = t27 * qJ(5);
t1 = -t3 + t80 + t82;
t94 = t27 * pkin(4);
t2 = -t4 - t94;
t6 = qJ(5) * t31 + t99;
t61 = t20 * t50 - t23 * t48;
t7 = -pkin(4) * t31 - t61;
t97 = t1 * t48 - t2 * t50 + (t48 * t7 + t50 * t6) * qJD(4);
t96 = 0.2e1 * qJD(4);
t95 = 0.2e1 * qJD(5);
t12 = t46 * t25 - t47 * t56;
t93 = t12 * t48;
t38 = pkin(2) * t46 + pkin(7);
t92 = t27 * t38;
t91 = t31 * t38;
t90 = t32 * t50;
t89 = t48 * t27;
t88 = t48 * t28;
t87 = t50 * t27;
t86 = t50 * t28;
t44 = t48 ^ 2;
t45 = t50 ^ 2;
t83 = t44 - t45;
t79 = t48 * qJD(5);
t75 = -0.2e1 * pkin(1) * qJD(2);
t39 = -pkin(2) * t47 - pkin(3);
t74 = t39 * t96;
t73 = t38 * t81;
t72 = t38 * t43;
t71 = t48 * t43;
t70 = -0.4e1 * t48 * t90;
t22 = -t47 * t34 - t46 * t35;
t67 = t83 * qJD(4);
t64 = -t48 * t6 + t50 * t7;
t62 = pkin(4) * t48 - qJ(5) * t50;
t59 = t28 * t39 - t92;
t58 = -t32 * t39 + t91;
t19 = t31 * t43 + t89;
t17 = -t31 * t81 + t87;
t57 = t32 * t43 + t88;
t18 = -t32 * t81 + t86;
t29 = t39 - t63;
t5 = t62 * t28 + t32 * t98 + t12;
t55 = -t5 + (t29 * t32 - t91) * qJD(4);
t26 = -pkin(4) * t81 + qJ(5) * t43 + t79;
t8 = t32 * t62 + t22;
t54 = qJD(4) * t8 - t26 * t32 + t28 * t29 - t92;
t52 = qJD(4) * t64 + t1 * t50 + t2 * t48;
t30 = t32 ^ 2;
t9 = [0, 0, 0, 0.2e1 * t49 * t76, 0.2e1 * (-t49 ^ 2 + t51 ^ 2) * qJD(2), 0, 0, 0, t49 * t75, t51 * t75, 0.2e1 * t27 * t40 + 0.2e1 * t31 * t42, 0.2e1 * t28 * t40 + 0.2e1 * t32 * t42, 0.2e1 * t12 * t32 - 0.2e1 * t13 * t31 + 0.2e1 * t22 * t28 - 0.2e1 * t23 * t27, 0.2e1 * t12 * t22 + 0.2e1 * t13 * t23 + 0.2e1 * t40 * t42, 0.2e1 * t28 * t32 * t45 - 0.2e1 * t30 * t71, t30 * t83 * t96 + t28 * t70, 0.2e1 * t18 * t31 + 0.2e1 * t32 * t87, -0.2e1 * t31 * t57 - 0.2e1 * t32 * t89, 0.2e1 * t31 * t27, 0.2e1 * t22 * t57 + 0.2e1 * t27 * t61 + 0.2e1 * t4 * t31 + 0.2e1 * t32 * t93, 0.2e1 * t12 * t90 + 0.2e1 * t18 * t22 - 0.2e1 * t27 * t99 + 0.2e1 * t3 * t31, 0.2e1 * t8 * t88 - 0.2e1 * t2 * t31 - 0.2e1 * t7 * t27 + 0.2e1 * (t43 * t8 + t5 * t48) * t32, 0.2e1 * t64 * t28 - 0.2e1 * t32 * t97, -0.2e1 * t8 * t86 + 0.2e1 * t1 * t31 + 0.2e1 * t6 * t27 + 0.2e1 * (-t5 * t50 + t8 * t81) * t32, 0.2e1 * t1 * t6 + 0.2e1 * t2 * t7 + 0.2e1 * t5 * t8; 0, 0, 0, 0, 0, t76, -t78, 0, -pkin(6) * t76, pkin(6) * t78, -t12, -t13, (-t27 * t46 - t28 * t47) * pkin(2), (-t12 * t47 + t13 * t46) * pkin(2), -t32 * t67 + t48 * t86, qJD(4) * t70 - t28 * t83, t19, t17, 0, -t12 * t50 + t59 * t48 + (t22 * t48 - t50 * t58) * qJD(4), t93 + t59 * t50 + (t22 * t50 + t48 * t58) * qJD(4), t48 * t54 + t50 * t55, t52, t48 * t55 - t50 * t54, -t8 * t26 + t5 * t29 + t38 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t71, -0.2e1 * t67, 0, 0, 0, t48 * t74, t50 * t74, 0.2e1 * t26 * t50 + 0.2e1 * t29 * t81, 0, 0.2e1 * t26 * t48 - 0.2e1 * t29 * t43, -0.2e1 * t29 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t28, 0, t42, 0, 0, 0, 0, 0, t17, -t19, t17, (-t44 - t45) * t28, t19, t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t57, t27, t4, t3, t4 + 0.2e1 * t94, -t63 * t28 + (qJD(4) * t62 - t79) * t32, -t3 + 0.2e1 * t80 + 0.2e1 * t82, -pkin(4) * t2 + qJ(5) * t1 + qJD(5) * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t81, 0, -t72, t73, -t72, -t98, -t73, -t98 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, -t43, -t81, 0, t43, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, qJ(5) * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, t18, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, 0, t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;
