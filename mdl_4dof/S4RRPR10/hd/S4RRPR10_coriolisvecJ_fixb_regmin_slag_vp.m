% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% tauc_reg [4x21]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRPR10_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR10_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:11:55
% EndTime: 2019-12-31 17:11:58
% DurationCPUTime: 0.69s
% Computational Cost: add. (434->145), mult. (1118->228), div. (0->0), fcn. (610->4), ass. (0->97)
t48 = sin(qJ(4));
t50 = cos(qJ(4));
t81 = t50 * qJD(2);
t51 = cos(qJ(2));
t89 = qJD(1) * t51;
t21 = -t48 * t89 + t81;
t78 = qJD(1) * qJD(2);
t107 = -0.2e1 * t78;
t104 = pkin(3) + pkin(5);
t49 = sin(qJ(2));
t106 = t104 * t49;
t83 = t49 * qJD(1);
t40 = qJD(4) + t83;
t105 = -qJD(4) + t40;
t52 = -pkin(2) - pkin(6);
t84 = t48 * qJD(2);
t19 = t50 * t89 + t84;
t69 = t49 * t78;
t4 = -qJD(4) * t19 + t48 * t69;
t103 = t4 * t50;
t11 = (-qJD(1) * t106 + qJD(3)) * qJD(2);
t102 = t11 * t48;
t101 = t11 * t50;
t100 = t19 * t40;
t99 = t21 * t40;
t98 = t21 * t51;
t97 = t40 * t48;
t96 = t40 * t49;
t95 = t40 * t52;
t54 = qJD(1) ^ 2;
t94 = t51 * t54;
t53 = qJD(2) ^ 2;
t93 = t53 * t49;
t92 = t53 * t51;
t46 = t49 ^ 2;
t47 = t51 ^ 2;
t91 = t46 - t47;
t90 = qJD(2) * pkin(2);
t67 = -t49 * qJ(3) - pkin(1);
t29 = -t51 * pkin(2) + t67;
t16 = qJD(1) * t29;
t88 = qJD(2) * t49;
t87 = qJD(2) * t51;
t86 = qJD(4) * t50;
t42 = pkin(5) * t89;
t79 = qJD(2) * qJ(3);
t30 = -t42 - t79;
t43 = pkin(3) * t89;
t15 = -t30 + t43;
t85 = t15 * qJD(4);
t82 = t49 * qJD(3);
t41 = pkin(5) * t83;
t80 = pkin(3) * t83 + qJD(3) + t41;
t77 = t50 * t96;
t76 = t49 * t94;
t32 = t104 * t51;
t74 = qJD(4) * t97;
t73 = t40 * t86;
t72 = t51 * t86;
t68 = t51 * t78;
t38 = pkin(5) * t68;
t17 = pkin(3) * t68 + t38;
t39 = pkin(2) * t69;
t64 = pkin(6) * t49 - qJ(3) * t51;
t56 = qJD(2) * t64 - t82;
t3 = qJD(1) * t56 + t39;
t70 = t50 * t17 - t48 * t3;
t66 = pkin(1) * t107;
t65 = qJD(3) - t90;
t10 = t52 * qJD(2) + t80;
t18 = t52 * t51 + t67;
t8 = t18 * qJD(1);
t1 = t50 * t10 - t48 * t8;
t2 = t48 * t10 + t50 * t8;
t63 = t106 * t48 + t50 * t18;
t62 = -qJD(1) * t47 + t96;
t61 = -0.2e1 * qJD(2) * t16;
t44 = pkin(2) * t88;
t57 = -t51 * t79 - t82;
t14 = t44 + t57;
t7 = qJD(1) * t57 + t39;
t59 = pkin(5) * t53 + qJD(1) * t14 + t7;
t5 = t21 * qJD(4) - t50 * t69;
t58 = t15 * t49 + t52 * t87;
t27 = (-qJD(3) + t41) * qJD(2);
t28 = t41 + t65;
t55 = -t27 * t51 + (t28 * t51 + (t30 + t42) * t49) * qJD(2);
t45 = pkin(2) * t83;
t34 = t50 * t68;
t26 = qJD(2) * t32;
t25 = t42 + t43;
t24 = qJD(2) * t106;
t22 = -qJ(3) * t89 + t45;
t13 = qJD(1) * t64 + t45;
t9 = t16 * t83;
t6 = t44 + t56;
t12 = [0, 0, 0, 0.2e1 * t49 * t68, t91 * t107, t92, -t93, 0, -pkin(5) * t92 + t49 * t66, pkin(5) * t93 + t51 * t66, t55, t49 * t61 + t51 * t59, -t49 * t59 + t51 * t61, pkin(5) * t55 + t16 * t14 + t7 * t29, -t4 * t48 * t51 + (t49 * t84 - t72) * t21, (-t19 * t48 + t21 * t50) * t88 + (-t103 + t48 * t5 + (t19 * t50 + t21 * t48) * qJD(4)) * t51, -t40 * t72 + t4 * t49 + (t48 * t62 + t98) * qJD(2), t51 * t74 - t5 * t49 + (-t19 * t51 + t50 * t62) * qJD(2), (t40 + t83) * t87, (t50 * t26 - t48 * t6) * t40 - t24 * t19 + t32 * t5 + (-t15 * t81 + t70) * t49 + (-t2 * t49 - t40 * t63) * qJD(4) + (-t48 * t85 + t101 + ((t106 * t50 - t48 * t18) * qJD(1) + t1) * qJD(2)) * t51, -t24 * t21 + t32 * t4 + (-(qJD(4) * t106 + t6) * t40 - (qJD(4) * t10 + t3) * t49) * t50 + (-(-qJD(4) * t18 + t26) * t40 + (t15 * qJD(2) + qJD(4) * t8 - t17) * t49) * t48 + (-t50 * t85 - t102 + (-qJD(1) * t63 - t2) * qJD(2)) * t51; 0, 0, 0, -t76, t91 * t54, 0, 0, 0, t54 * pkin(1) * t49, pkin(1) * t94, ((-t30 - t79) * t49 + (-t28 + t65) * t51) * qJD(1), -t22 * t89 + t9, 0.2e1 * qJD(2) * qJD(3) + (t16 * t51 + t22 * t49) * qJD(1), -t27 * qJ(3) - t30 * qJD(3) - t16 * t22 + (-t30 * t49 + (-t28 - t90) * t51) * qJD(1) * pkin(5), -t21 * t97 + t103, (-t5 - t99) * t50 + (-t4 + t100) * t48, -t74 + t34 + (-t48 * t96 - t98) * qJD(1), -t73 + (-t77 + (t19 - t84) * t51) * qJD(1), -t40 * t89, qJ(3) * t5 + t102 - (-t48 * t13 + t50 * t25) * t40 + t80 * t19 + (t15 * t50 - t48 * t95) * qJD(4) + (-t1 * t51 + t50 * t58) * qJD(1), qJ(3) * t4 + t101 + (t50 * t13 + t48 * t25) * t40 + t80 * t21 + (-t15 * t48 - t50 * t95) * qJD(4) + (t2 * t51 - t48 * t58) * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, -t46 * t54 - t53, t30 * qJD(2) + t38 + t9, 0, 0, 0, 0, 0, -qJD(2) * t19 - t40 * t97 + t34, -t73 - qJD(2) * t21 + (-t51 * t84 - t77) * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21 * t19, -t19 ^ 2 + t21 ^ 2, t4 + t100, -t5 + t99, t68, t105 * t2 - t15 * t21 + t70, t105 * t1 + t15 * t19 - t48 * t17 - t50 * t3;];
tauc_reg = t12;
