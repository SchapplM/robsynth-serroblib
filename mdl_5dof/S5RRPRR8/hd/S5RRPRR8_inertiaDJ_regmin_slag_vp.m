% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPRR8
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
% Datum: 2019-12-31 20:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRR8_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR8_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:18:13
% EndTime: 2019-12-31 20:18:15
% DurationCPUTime: 0.66s
% Computational Cost: add. (1365->115), mult. (3126->220), div. (0->0), fcn. (3100->8), ass. (0->88)
t106 = cos(qJ(4));
t62 = sin(pkin(9));
t108 = t62 * pkin(2);
t65 = sin(qJ(4));
t63 = cos(pkin(9));
t87 = t63 * pkin(2) + pkin(3);
t109 = -t106 * t87 + t65 * t108;
t67 = cos(qJ(5));
t61 = t67 ^ 2;
t64 = sin(qJ(5));
t98 = t64 ^ 2 - t61;
t82 = t98 * qJD(5);
t100 = -qJ(3) - pkin(6);
t66 = sin(qJ(2));
t54 = t100 * t66;
t68 = cos(qJ(2));
t55 = t100 * t68;
t29 = t63 * t54 + t62 * t55;
t48 = t62 * t68 + t63 * t66;
t24 = -t48 * pkin(7) + t29;
t30 = t62 * t54 - t63 * t55;
t47 = -t62 * t66 + t63 * t68;
t25 = t47 * pkin(7) + t30;
t83 = qJD(2) * t100;
t43 = t68 * qJD(3) + t66 * t83;
t44 = -t66 * qJD(3) + t68 * t83;
t23 = t63 * t43 + t62 * t44;
t73 = t48 * qJD(2);
t69 = -pkin(7) * t73 + t23;
t22 = -t62 * t43 + t63 * t44;
t94 = t68 * qJD(2);
t95 = t66 * qJD(2);
t46 = -t62 * t95 + t63 * t94;
t72 = t46 * pkin(7) - t22;
t84 = qJD(4) * t106;
t97 = qJD(4) * t65;
t5 = t106 * t72 + t24 * t97 + t25 * t84 + t65 * t69;
t3 = t5 * t64;
t14 = -t106 * t24 + t65 * t25;
t59 = qJD(5) * t67;
t107 = t14 * t59 + t3;
t74 = t106 * t47 - t65 * t48;
t16 = t74 * qJD(4) + t106 * t46 - t65 * t73;
t28 = t106 * t48 + t65 * t47;
t105 = t28 * t16;
t104 = t28 * t67;
t17 = t28 * qJD(4) + t106 * t73 + t65 * t46;
t103 = t64 * t17;
t102 = t67 * t16;
t101 = t67 * t17;
t71 = t106 * t108 + t65 * t87;
t39 = t71 * qJD(4);
t41 = -pkin(4) + t109;
t99 = t39 * t64 + t41 * t59;
t96 = qJD(5) * t64;
t92 = -0.2e1 * pkin(1) * qJD(2);
t91 = pkin(4) * t96;
t90 = pkin(4) * t59;
t58 = pkin(2) * t95;
t89 = t64 * t59;
t88 = -t68 * pkin(2) - pkin(1);
t86 = -0.4e1 * t64 * t104;
t85 = -t39 * t67 + t41 * t96;
t35 = -t47 * pkin(3) + t88;
t13 = -pkin(4) * t74 - t28 * pkin(8) + t35;
t15 = t106 * t25 + t65 * t24;
t81 = t67 * t13 - t64 * t15;
t80 = t64 * t13 + t67 * t15;
t42 = pkin(8) + t71;
t79 = -t28 * t41 - t42 * t74;
t77 = t64 * t16 + t28 * t59;
t76 = t28 * t96 - t102;
t10 = -t59 * t74 + t103;
t75 = -t74 * t96 - t101;
t31 = pkin(3) * t73 + t58;
t38 = t109 * qJD(4);
t70 = t16 * t41 - t17 * t42 + t28 * t39 - t38 * t74;
t56 = 0.2e1 * t89;
t52 = -0.2e1 * t82;
t26 = t28 ^ 2;
t11 = t14 * t96;
t8 = t64 * t102 - t28 * t82;
t7 = t17 * pkin(4) - t16 * pkin(8) + t31;
t6 = qJD(5) * t86 - t98 * t16;
t4 = -t106 * t69 - t24 * t84 + t25 * t97 + t65 * t72;
t2 = -t80 * qJD(5) + t64 * t4 + t67 * t7;
t1 = -t81 * qJD(5) + t67 * t4 - t64 * t7;
t9 = [0, 0, 0, 0.2e1 * t66 * t94, 0.2e1 * (-t66 ^ 2 + t68 ^ 2) * qJD(2), 0, 0, 0, t66 * t92, t68 * t92, -0.2e1 * t22 * t48 + 0.2e1 * t23 * t47 - 0.2e1 * t29 * t46 - 0.2e1 * t30 * t73, 0.2e1 * t29 * t22 + 0.2e1 * t30 * t23 + 0.2e1 * t88 * t58, 0.2e1 * t105, 0.2e1 * t16 * t74 - 0.2e1 * t28 * t17, 0, 0, 0, 0.2e1 * t35 * t17 - 0.2e1 * t31 * t74, 0.2e1 * t35 * t16 + 0.2e1 * t31 * t28, 0.2e1 * t61 * t105 - 0.2e1 * t26 * t89, t16 * t86 + 0.2e1 * t26 * t82, 0.2e1 * t28 * t101 + 0.2e1 * t74 * t76, -0.2e1 * t28 * t103 + 0.2e1 * t74 * t77, -0.2e1 * t74 * t17, 0.2e1 * t77 * t14 + 0.2e1 * t81 * t17 - 0.2e1 * t2 * t74 + 0.2e1 * t28 * t3, -0.2e1 * t1 * t74 + 0.2e1 * t5 * t104 - 0.2e1 * t76 * t14 - 0.2e1 * t80 * t17; 0, 0, 0, 0, 0, t94, -t95, 0, -pkin(6) * t94, pkin(6) * t95, (-t63 * t46 - t62 * t73) * pkin(2), (t22 * t63 + t23 * t62) * pkin(2), 0, 0, t16, -t17, 0, -t5, t4, t8, t6, t10, -t75, 0, t11 + (-t79 * qJD(5) - t5) * t67 + t70 * t64, t70 * t67 + t79 * t96 + t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t39, 0.2e1 * t38, t56, t52, 0, 0, 0, 0.2e1 * t85, 0.2e1 * t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, 0, 0, 0, 0, t17, t16, 0, 0, 0, 0, 0, -t75, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t17, 0, -t5, t4, t8, t6, t10, -t75, 0, t11 + (-pkin(4) * t16 - pkin(8) * t17) * t64 + (-t5 + (-pkin(4) * t28 + pkin(8) * t74) * qJD(5)) * t67, pkin(4) * t76 + pkin(8) * t75 + t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, t38, t56, t52, 0, 0, 0, t85 - t91, -t90 + t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, t52, 0, 0, 0, -0.2e1 * t91, -0.2e1 * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, -t77, t17, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, -t96, 0, t64 * t38 - t42 * t59, t67 * t38 + t42 * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, -t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, -t96, 0, -pkin(8) * t59, pkin(8) * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;
