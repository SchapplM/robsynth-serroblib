% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRP4_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP4_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:52:56
% EndTime: 2019-12-31 19:52:58
% DurationCPUTime: 0.56s
% Computational Cost: add. (730->140), mult. (1179->160), div. (0->0), fcn. (480->4), ass. (0->90)
t52 = sin(qJ(4));
t50 = t52 ^ 2;
t54 = cos(qJ(4));
t51 = t54 ^ 2;
t90 = t50 + t51;
t110 = 2 * qJD(4);
t84 = qJD(4) * qJ(5);
t49 = qJD(1) + qJD(2);
t56 = -pkin(2) - pkin(7);
t55 = cos(qJ(2));
t89 = pkin(1) * qJD(1);
t80 = t55 * t89;
t68 = qJD(3) - t80;
t16 = t56 * t49 + t68;
t97 = t52 * t16;
t10 = t84 + t97;
t53 = sin(qJ(2));
t88 = pkin(1) * qJD(2);
t78 = qJD(1) * t88;
t71 = t53 * t78;
t36 = t52 * t71;
t96 = t54 * t16;
t3 = t36 + (qJD(5) + t96) * qJD(4);
t37 = t54 * t71;
t86 = qJD(4) * t52;
t5 = t16 * t86 - t37;
t74 = (qJD(4) * pkin(4)) - qJD(5);
t9 = -t74 - t96;
t58 = t3 * t52 - t5 * t54 + (t10 * t54 + t52 * t9) * qJD(4);
t35 = t52 * pkin(4) - t54 * qJ(5) + qJ(3);
t109 = t35 * t49;
t83 = t53 * t88;
t63 = t90 * t83;
t108 = qJD(1) + t49;
t48 = t49 ^ 2;
t91 = t50 - t51;
t107 = t91 * t49 * t110;
t41 = t55 * t78;
t65 = pkin(4) * t54 + qJ(5) * t52;
t72 = -t54 * qJD(5) + qJD(3);
t2 = t41 + (qJD(4) * t65 + t72) * t49;
t81 = t53 * t89;
t8 = t81 + t109;
t85 = qJD(4) * t54;
t106 = t2 * t52 + t8 * t85;
t104 = t53 * pkin(1);
t87 = t49 * qJ(3);
t30 = t81 + t87;
t103 = t30 * t49;
t102 = t30 * t55;
t79 = -t55 * pkin(1) - pkin(2);
t42 = -pkin(7) + t79;
t57 = qJD(4) ^ 2;
t101 = t42 * t57;
t100 = t48 * t54;
t99 = t49 * t52;
t98 = t49 * t54;
t95 = t56 * t57;
t94 = t57 * t52;
t26 = t49 * qJD(3) + t41;
t93 = t26 * t52 + t30 * t85;
t92 = t90 * t49 * t81;
t82 = t55 * t88;
t76 = t85 * t104;
t73 = t85 * t99;
t70 = qJD(2) * t76 - t42 * t94;
t15 = pkin(4) * t85 + t52 * t84 + t72;
t69 = -t15 + t80;
t66 = t10 * t52 - t54 * t9;
t40 = qJD(3) + t82;
t43 = qJ(3) + t104;
t64 = t26 * t43 + t30 * t40;
t62 = t26 * qJ(3) + t30 * qJD(3);
t61 = -qJD(1) * t76 - t56 * t94;
t46 = t57 * t54;
t39 = t52 * t100;
t34 = -0.2e1 * t73;
t33 = 0.2e1 * t73;
t29 = t46 + t100;
t28 = (-t48 - t57) * t52;
t27 = -t49 * pkin(2) + t68;
t25 = t91 * t48;
t24 = t35 + t104;
t22 = t65 * t49;
t21 = t26 * t54;
t19 = t108 * t83;
t18 = (qJD(2) - t49) * t81;
t11 = t15 + t82;
t6 = t8 * t86;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t49 * t82 - t41, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, t41 + (qJD(3) + t40) * t49, (t79 * qJD(1) + t27) * t83 + t64, t34, t107, -t94, t33, -t46, 0, (t40 * t52 + t43 * t85) * t49 + t70 + t93, t21 + (t40 * t49 - t101) * t54 + (-t43 * t49 - t30 - t83) * t86, -t108 * t63, t64 + (qJD(1) * t42 + t16) * t63, t34, -t94, -t107, 0, t46, t33, (t11 * t52 + t24 * t85) * t49 + t70 + t106, -t49 * t63 - t58, t6 + (t24 * t49 + t83) * t86 + (-t11 * t49 + t101 - t2) * t54, t8 * t11 + t2 * t24 + t42 * t58 + t66 * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, t49 * t80 - t41, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, t41 + (0.2e1 * qJD(3) - t80) * t49, (-t102 + (-pkin(2) * qJD(2) - t27) * t53) * t89 + t62, t34, t107, -t94, t33, -t46, 0, (qJ(3) * t85 + t52 * t68) * t49 + t61 + t93, t21 + (t49 * t68 - t95) * t54 + (-t30 + t81 - t87) * t86, -qJD(1) * t63 + t92, (-t102 + (qJD(2) * t56 - t16) * t53 * t90) * t89 + t62, t34, -t94, -t107, 0, t46, t33, (t35 * t85 - t52 * t69) * t49 + t61 + t106, -t58 + t92, t6 + (-t81 + t109) * t86 + (t49 * t69 - t2 + t95) * t54, t8 * t15 + t2 * t35 + (-t53 * t66 - t55 * t8) * t89 + t58 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, t71 - t103, 0, 0, 0, 0, 0, 0, t28, -t29, 0, t90 * t71 - t103, 0, 0, 0, 0, 0, 0, t28, 0, t29, -t8 * t49 + t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t25, 0, -t39, 0, 0, -t30 * t98 + t37, t30 * t99 - t36, 0, 0, t39, 0, t25, 0, 0, -t39, t37 + (-t22 * t52 - t54 * t8) * t49, ((t10 - t84) * t54 + (t74 + t9) * t52) * t49, qJD(5) * t110 + t36 + (t22 * t54 - t52 * t8) * t49, -t9 * t97 - t5 * pkin(4) + t3 * qJ(5) - t8 * t22 + (qJD(5) - t96) * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, -t51 * t48 - t57, -t10 * qJD(4) + t8 * t98 + t5;];
tauc_reg = t1;
