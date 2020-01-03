% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
% 
% Output:
% tauc_reg [5x21]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPP5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:42
% EndTime: 2019-12-31 18:16:44
% DurationCPUTime: 0.45s
% Computational Cost: add. (455->111), mult. (882->151), div. (0->0), fcn. (344->2), ass. (0->86)
t47 = -pkin(3) - pkin(4);
t48 = -pkin(1) - pkin(6);
t32 = t48 * qJD(1) + qJD(2);
t46 = cos(qJ(3));
t74 = qJ(5) * qJD(1);
t16 = (t32 + t74) * t46;
t77 = qJD(4) - t16;
t7 = qJD(3) * t47 + t77;
t92 = t46 * t32;
t24 = qJD(3) * t92;
t41 = qJD(3) * qJD(4);
t18 = t24 + t41;
t85 = qJD(3) * pkin(3);
t61 = -qJD(4) + t85;
t19 = -t61 - t92;
t45 = sin(qJ(3));
t26 = t45 * t32;
t42 = qJD(3) * qJ(4);
t21 = t26 + t42;
t93 = t21 * t46;
t98 = ((-t19 + t92) * t45 - t93) * qJD(3) - t18 * t45;
t70 = 2 * qJD(1);
t64 = t46 * qJ(4) - qJ(2);
t22 = t45 * t47 + t64;
t82 = qJD(1) * t22;
t10 = qJD(5) + t82;
t97 = (qJD(5) + t10) * t46;
t49 = qJD(3) ^ 2;
t91 = t49 * t45;
t90 = t49 * t46;
t81 = qJD(3) * t45;
t23 = t32 * t81;
t73 = qJD(1) * qJD(3);
t65 = qJ(5) * t73;
t89 = t45 * t65 + t23;
t43 = t45 ^ 2;
t44 = t46 ^ 2;
t88 = t43 - t44;
t50 = qJD(1) ^ 2;
t87 = t49 + t50;
t86 = qJ(4) * t45;
t84 = t50 * qJ(2);
t83 = qJ(5) + t48;
t29 = t45 * pkin(3) - t64;
t20 = qJD(1) * t29;
t39 = t45 * t74;
t11 = t39 + t21;
t80 = t11 * qJD(3);
t79 = t20 * qJD(1);
t78 = t46 * qJD(4);
t75 = qJ(2) * qJD(3);
t72 = qJD(1) * qJD(5);
t71 = t45 * t72 + t46 * t65 + t24;
t69 = qJD(2) * t70;
t68 = -0.2e1 * t73;
t38 = qJD(1) * t78;
t54 = t47 * t46 - t86;
t52 = t54 * qJD(3) - qJD(2);
t2 = t52 * qJD(1) + t38;
t9 = t52 + t78;
t67 = qJD(1) * t9 + t2;
t63 = -t10 - t82;
t62 = qJD(3) * t83;
t59 = t46 * t68;
t3 = t41 + t71;
t5 = -t46 * t72 + t89;
t58 = t3 * t45 + t7 * t81 + (-t5 + t80) * t46;
t57 = pkin(3) * t46 + t86;
t56 = 0.2e1 * qJD(3) * t20;
t53 = t57 * qJD(3) + qJD(2);
t14 = t53 - t78;
t8 = t53 * qJD(1) - t38;
t55 = qJD(1) * t14 - t48 * t49 + t8;
t40 = 0.2e1 * t41;
t36 = t46 * t50 * t45;
t33 = -t44 * t50 - t49;
t31 = t87 * t46;
t30 = t87 * t45;
t28 = t83 * t46;
t27 = t83 * t45;
t25 = t57 * qJD(1);
t17 = t54 * qJD(1);
t15 = t26 + t39;
t13 = t45 * qJD(5) + t46 * t62;
t12 = -t46 * qJD(5) + t45 * t62;
t1 = [0, 0, 0, 0, t69, qJ(2) * t69, t45 * t59, 0.2e1 * t88 * t73, -t91, -t90, 0, -t48 * t91 + (qJD(2) * t45 + t46 * t75) * t70, -t48 * t90 + (qJD(2) * t46 - t45 * t75) * t70, t55 * t45 + t46 * t56, t98, t45 * t56 - t55 * t46, t20 * t14 + t8 * t29 - t48 * t98, -t67 * t45 + (t63 * t46 - t12) * qJD(3), t67 * t46 + (t63 * t45 + t13) * qJD(3), (-t12 * t46 + t13 * t45 + (t27 * t46 - t28 * t45) * qJD(3)) * qJD(1) + t58, t10 * t9 + t11 * t13 + t7 * t12 + t2 * t22 + t3 * t27 - t5 * t28; 0, 0, 0, 0, -t50, -t84, 0, 0, 0, 0, 0, -t30, -t31, -t30, 0, t31, -t98 - t79, -t30, t31, 0, t10 * qJD(1) + t58; 0, 0, 0, 0, 0, 0, t36, -t88 * t50, 0, 0, 0, -t46 * t84, t45 * t84, (-t20 * t46 - t25 * t45) * qJD(1), ((t21 - t42) * t46 + (t19 + t61) * t45) * qJD(1), t40 + (-t20 * t45 + t25 * t46) * qJD(1), t18 * qJ(4) + t21 * qJD(4) - t20 * t25 + (-t93 + (-t19 - t85) * t45) * t32, t15 * qJD(3) + (t17 * t45 + t97) * qJD(1) - t89, -t16 * qJD(3) + t40 + (t10 * t45 - t17 * t46) * qJD(1) + t71, (-t11 + t15 + t42) * t46 * qJD(1), t3 * qJ(4) - t10 * t17 + t77 * t11 - t7 * t15 + t5 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, t33, -t21 * qJD(3) + t46 * t79 + t23, t36, t33, 0, -qJD(1) * t97 - t80 + t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, t45 * t68, (-t43 - t44) * t50, t38 + (-t11 * t45 + t46 * t7 + t52) * qJD(1);];
tauc_reg = t1;
