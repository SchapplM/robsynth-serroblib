% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% tauc_reg [4x19]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRPR7_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR7_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:06:44
% EndTime: 2019-12-31 17:06:46
% DurationCPUTime: 0.64s
% Computational Cost: add. (804->144), mult. (2203->222), div. (0->0), fcn. (1535->6), ass. (0->92)
t63 = cos(pkin(7));
t67 = cos(qJ(2));
t102 = t63 * t67;
t53 = qJD(1) * t102;
t62 = sin(pkin(7));
t65 = sin(qJ(2));
t91 = t65 * qJD(1);
t37 = -t62 * t91 + t53;
t35 = qJD(4) - t37;
t113 = qJD(4) - t35;
t89 = qJD(1) * qJD(2);
t112 = -0.2e1 * t89;
t46 = t62 * t67 + t63 * t65;
t39 = t46 * qJD(1);
t96 = -qJ(3) - pkin(5);
t83 = qJD(2) * t96;
t36 = t67 * qJD(3) + t65 * t83;
t30 = t36 * qJD(1);
t72 = -t65 * qJD(3) + t67 * t83;
t70 = t72 * qJD(1);
t4 = t30 * t62 - t63 * t70;
t56 = pkin(2) * t62 + pkin(6);
t111 = (pkin(2) * t91 + pkin(3) * t39 - pkin(6) * t37 + qJD(4) * t56) * t35 + t4;
t64 = sin(qJ(4));
t66 = cos(qJ(4));
t26 = qJD(2) * t64 + t39 * t66;
t85 = t65 * t89;
t33 = qJD(2) * t53 - t62 * t85;
t7 = qJD(4) * t26 + t64 * t33;
t11 = t63 * t36 + t62 * t72;
t51 = t96 * t67;
t49 = qJD(1) * t51;
t103 = t62 * t49;
t86 = t96 * t65;
t48 = qJD(1) * t86;
t94 = qJD(2) * pkin(2);
t44 = t48 + t94;
t18 = t44 * t63 + t103;
t14 = -qJD(2) * pkin(3) - t18;
t45 = t62 * t65 - t102;
t87 = -pkin(2) * t67 - pkin(1);
t17 = pkin(3) * t45 - pkin(6) * t46 + t87;
t41 = t45 * qJD(2);
t5 = t63 * t30 + t62 * t70;
t23 = -t63 * t51 + t62 * t86;
t38 = t46 * qJD(2);
t32 = qJD(1) * t38;
t77 = -t23 * t32 + t4 * t46;
t76 = t87 * qJD(1);
t50 = qJD(3) + t76;
t9 = -t37 * pkin(3) - t39 * pkin(6) + t50;
t110 = -t14 * t41 - (qJD(4) * t17 + t11) * t35 - (qJD(4) * t9 + t5) * t45 + t77;
t90 = t66 * qJD(2);
t92 = qJD(4) * t64;
t6 = qJD(4) * t90 + t33 * t66 - t39 * t92;
t109 = t6 * t64;
t108 = t17 * t32;
t24 = t39 * t64 - t90;
t107 = t24 * t35;
t106 = t26 * t35;
t105 = t26 * t39;
t104 = t39 * t24;
t42 = t63 * t49;
t101 = t64 * t32;
t28 = t66 * t32;
t69 = qJD(1) ^ 2;
t99 = t67 * t69;
t68 = qJD(2) ^ 2;
t98 = t68 * t65;
t97 = t68 * t67;
t19 = t44 * t62 - t42;
t95 = t65 ^ 2 - t67 ^ 2;
t93 = qJD(4) * t46;
t88 = t65 * t94;
t82 = t35 * t66;
t79 = pkin(1) * t112;
t15 = qJD(2) * pkin(6) + t19;
t2 = t15 * t66 + t64 * t9;
t78 = t15 * t64 - t66 * t9;
t75 = t28 + (t37 * t64 - t92) * t35;
t74 = -t41 * t66 - t46 * t92;
t21 = t48 * t63 + t103;
t71 = -t56 * t32 + (t14 + t21) * t35;
t57 = -pkin(2) * t63 - pkin(3);
t54 = pkin(2) * t85;
t22 = -t51 * t62 - t63 * t86;
t20 = t48 * t62 - t42;
t13 = pkin(3) * t38 + pkin(6) * t41 + t88;
t10 = t36 * t62 - t63 * t72;
t8 = pkin(3) * t32 - pkin(6) * t33 + t54;
t3 = t66 * t8;
t1 = [0, 0, 0, 0.2e1 * t67 * t85, t95 * t112, t97, -t98, 0, -pkin(5) * t97 + t65 * t79, pkin(5) * t98 + t67 * t79, t10 * t39 + t11 * t37 + t18 * t41 - t19 * t38 + t22 * t33 - t5 * t45 + t77, -t18 * t10 + t19 * t11 + t4 * t22 + t5 * t23 + (t50 + t76) * t88, t6 * t66 * t46 + t26 * t74, -(-t24 * t66 - t26 * t64) * t41 + (-t109 - t66 * t7 + (t24 * t64 - t26 * t66) * qJD(4)) * t46, t26 * t38 + t28 * t46 + t35 * t74 + t6 * t45, -t46 * t101 - t24 * t38 - t7 * t45 + (t41 * t64 - t66 * t93) * t35, t32 * t45 + t35 * t38, -t78 * t38 + t10 * t24 + t22 * t7 + t3 * t45 + (t13 * t35 + t108 + (t14 * t46 - t15 * t45 - t23 * t35) * qJD(4)) * t66 + t110 * t64, t10 * t26 - t2 * t38 + t22 * t6 + (-(-qJD(4) * t23 + t13) * t35 - t108 - (-qJD(4) * t15 + t8) * t45 - t14 * t93) * t64 + t110 * t66; 0, 0, 0, -t65 * t99, t95 * t69, 0, 0, 0, t69 * pkin(1) * t65, pkin(1) * t99, (t19 - t20) * t39 + (t18 - t21) * t37 + (-t32 * t62 - t33 * t63) * pkin(2), t18 * t20 - t19 * t21 + (-t4 * t63 + t5 * t62 - t50 * t91) * pkin(2), t26 * t82 + t109, (t6 - t107) * t66 + (-t7 - t106) * t64, t35 * t82 + t101 - t105, t75 + t104, -t35 * t39, -t111 * t66 - t20 * t24 + t39 * t78 + t57 * t7 + t71 * t64, t111 * t64 + t2 * t39 - t20 * t26 + t57 * t6 + t71 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37 ^ 2 - t39 ^ 2, t18 * t39 - t19 * t37 + t54, 0, 0, 0, 0, 0, t75 - t104, -t35 ^ 2 * t66 - t101 - t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26 * t24, -t24 ^ 2 + t26 ^ 2, t6 + t107, t106 - t7, t32, -t113 * t2 - t14 * t26 - t64 * t5 + t3, t113 * t78 + t14 * t24 - t66 * t5 - t64 * t8;];
tauc_reg = t1;
