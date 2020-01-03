% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% tauc_reg [4x21]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPRR6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:52:41
% EndTime: 2019-12-31 16:52:43
% DurationCPUTime: 0.48s
% Computational Cost: add. (614->105), mult. (1747->160), div. (0->0), fcn. (1354->6), ass. (0->72)
t54 = cos(pkin(7));
t58 = cos(qJ(3));
t83 = t58 * t54;
t53 = sin(pkin(7));
t56 = sin(qJ(3));
t85 = t56 * t53;
t67 = -t83 + t85;
t35 = t67 * qJD(1);
t57 = cos(qJ(4));
t29 = t57 * t35;
t79 = t58 * qJD(3);
t47 = t54 * qJD(1) * t79;
t77 = qJD(1) * t85;
t32 = -qJD(3) * t77 + t47;
t41 = t58 * t53 + t56 * t54;
t38 = t41 * qJD(3);
t33 = qJD(1) * t38;
t36 = t41 * qJD(1);
t55 = sin(qJ(4));
t80 = qJD(4) * t55;
t1 = -qJD(4) * t29 + t57 * t32 - t55 * t33 - t36 * t80;
t17 = t55 * t36 + t29;
t52 = qJD(3) + qJD(4);
t87 = t17 * t52;
t98 = t1 + t87;
t71 = -t55 * t35 + t57 * t36;
t97 = t71 * t17;
t2 = t71 * qJD(4) + t55 * t32 + t57 * t33;
t88 = t71 * t52;
t96 = -t2 + t88;
t95 = -t17 ^ 2 + t71 ^ 2;
t82 = pkin(5) + qJ(2);
t45 = t82 * t53;
t42 = qJD(1) * t45;
t46 = t82 * t54;
t43 = qJD(1) * t46;
t69 = t56 * t42 - t58 * t43;
t13 = -t35 * pkin(6) - t69;
t49 = -t54 * pkin(2) - pkin(1);
t44 = t49 * qJD(1) + qJD(2);
t23 = t35 * pkin(3) + t44;
t64 = t41 * qJD(2);
t63 = qJD(1) * t64;
t7 = -t32 * pkin(6) + t69 * qJD(3) - t63;
t94 = t23 * t17 + t13 * t80 + (-t13 * t52 - t7) * t55;
t92 = -t58 * t42 - t56 * t43;
t91 = qJD(4) - t52;
t6 = -t33 * pkin(6) - qJD(2) * t35 + t92 * qJD(3);
t90 = -t23 * t71 - t55 * t6 + t57 * t7;
t89 = pkin(3) * t36;
t84 = t57 * t13;
t81 = t53 ^ 2 + t54 ^ 2;
t78 = qJD(1) * qJD(2);
t12 = -t36 * pkin(6) + t92;
t11 = qJD(3) * pkin(3) + t12;
t75 = -pkin(3) * t52 - t11;
t74 = t81 * qJD(1) ^ 2;
t70 = -t55 * t41 - t57 * t67;
t22 = t57 * t41 - t55 * t67;
t68 = t56 * t45 - t58 * t46;
t66 = 0.2e1 * t81 * t78;
t62 = -t45 * t79 + qJD(2) * t83 + (-qJD(2) * t53 - qJD(3) * t46) * t56;
t60 = t68 * qJD(3) - t64;
t37 = t67 * qJD(3);
t27 = pkin(3) * t67 + t49;
t15 = -pkin(6) * t67 - t68;
t14 = -t41 * pkin(6) - t58 * t45 - t56 * t46;
t9 = t37 * pkin(6) + t60;
t8 = -t38 * pkin(6) + t62;
t5 = t22 * qJD(4) - t55 * t37 + t57 * t38;
t4 = t70 * qJD(4) - t57 * t37 - t55 * t38;
t3 = [0, 0, 0, 0, 0, t66, qJ(2) * t66, t32 * t41 - t36 * t37, -t32 * t67 - t41 * t33 + t37 * t35 - t36 * t38, -t37 * qJD(3), -t38 * qJD(3), 0, t60 * qJD(3) + t49 * t33 + t44 * t38, -t62 * qJD(3) + t49 * t32 - t44 * t37, t1 * t22 + t4 * t71, t1 * t70 - t4 * t17 - t22 * t2 - t5 * t71, t4 * t52, -t5 * t52, 0, t27 * t2 + t23 * t5 + (-t55 * t8 + t57 * t9 + (-t14 * t55 - t15 * t57) * qJD(4)) * t52 + (t38 * t17 - t33 * t70) * pkin(3), t27 * t1 + t23 * t4 - (t55 * t9 + t57 * t8 + (t14 * t57 - t15 * t55) * qJD(4)) * t52 + (t33 * t22 + t38 * t71) * pkin(3); 0, 0, 0, 0, 0, -t74, -qJ(2) * t74, 0, 0, 0, 0, 0, 0.2e1 * t36 * qJD(3), t47 + (-t35 - t77) * qJD(3), 0, 0, 0, 0, 0, t2 + t88, t1 - t87; 0, 0, 0, 0, 0, 0, 0, t36 * t35, -t35 ^ 2 + t36 ^ 2, t47 + (t35 - t77) * qJD(3), 0, 0, -t44 * t36 - t63, t44 * t35 + t67 * t78, t97, t95, t98, t96, 0, -t17 * t89 - (-t55 * t12 - t84) * t52 + (t75 * t55 - t84) * qJD(4) + t90, -t71 * t89 + (t75 * qJD(4) + t12 * t52 - t6) * t57 + t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, t95, t98, t96, 0, t91 * (-t55 * t11 - t84) + t90, (-t91 * t11 - t6) * t57 + t94;];
tauc_reg = t3;
