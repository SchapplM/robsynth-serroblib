% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% 
% Output:
% tauc_reg [5x21]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPRP6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:55:16
% EndTime: 2019-12-31 17:55:18
% DurationCPUTime: 0.40s
% Computational Cost: add. (701->121), mult. (1534->152), div. (0->0), fcn. (987->4), ass. (0->77)
t56 = cos(pkin(7));
t59 = cos(qJ(4));
t85 = t59 * t56;
t74 = qJD(1) * t85;
t58 = sin(qJ(4));
t55 = sin(pkin(7));
t80 = qJD(1) * t55;
t75 = t58 * t80;
t29 = t74 - t75;
t78 = qJD(4) * t59;
t79 = qJD(4) * t58;
t31 = -t55 * t78 - t56 * t79;
t25 = t31 * qJD(4);
t34 = t59 * t55 + t58 * t56;
t64 = qJD(1) * t34;
t94 = -qJD(1) * t64 + t25;
t93 = 0.2e1 * t64;
t57 = -pkin(1) - qJ(3);
t92 = qJD(1) * t57;
t84 = t55 ^ 2 + t56 ^ 2;
t91 = t84 * qJD(3);
t90 = t29 ^ 2;
t53 = qJD(1) * qJD(2);
t73 = 0.2e1 * t53;
t54 = qJD(1) * qJ(2);
t49 = qJD(3) + t54;
t37 = pkin(3) * t80 + t49;
t9 = pkin(4) * t64 - t29 * qJ(5) + t37;
t89 = t9 * t29;
t88 = -pkin(6) + t57;
t87 = t29 * t64;
t41 = qJD(2) + t92;
t71 = -pkin(6) * qJD(1) + t41;
t23 = t71 * t55;
t86 = t58 * t23;
t35 = t88 * t55;
t36 = t88 * t56;
t16 = t58 * t35 - t59 * t36;
t65 = t34 * qJD(3);
t6 = -t16 * qJD(4) - t65;
t83 = t6 * qJD(4);
t17 = t59 * t35 + t58 * t36;
t33 = t58 * t55 - t85;
t7 = -t33 * qJD(3) + t17 * qJD(4);
t82 = t7 * qJD(4);
t46 = t55 * pkin(3) + qJ(2);
t32 = -t55 * t79 + t56 * t78;
t77 = t32 * qJD(4);
t24 = t71 * t56;
t10 = t59 * t24 - t86;
t76 = qJD(5) - t10;
t72 = t10 + t86;
t70 = qJD(1) * t84;
t2 = t29 * qJD(3) + t23 * t78 + t24 * t79;
t21 = qJD(4) * t64;
t69 = t33 * t21 + t29 * t31;
t11 = t59 * t23 + t58 * t24;
t68 = qJD(1) * t29 + t77;
t40 = qJD(4) * t75;
t22 = qJD(4) * t74 - t40;
t67 = t22 * pkin(4) + t21 * qJ(5) + t53;
t66 = t11 * qJD(4) - t2;
t63 = -t40 + (t29 + t74) * qJD(4);
t20 = t24 * t78;
t61 = -qJD(1) * t65 + t20;
t1 = (qJD(5) - t86) * qJD(4) + t61;
t5 = -qJD(4) * pkin(4) + t76;
t8 = qJD(4) * qJ(5) + t11;
t62 = t1 * t34 + t2 * t33 - t5 * t31 + t8 * t32;
t60 = qJD(1) ^ 2;
t26 = t64 ^ 2;
t15 = t29 * pkin(4) + qJ(5) * t64;
t14 = t34 * pkin(4) + t33 * qJ(5) + t46;
t12 = t93 * qJD(4);
t4 = t32 * pkin(4) - t31 * qJ(5) + t33 * qJD(5) + qJD(2);
t3 = -t29 * qJD(5) + t67;
t13 = [0, 0, 0, 0, t73, qJ(2) * t73, t55 * t73, t56 * t73, 0.2e1 * qJD(3) * t70, (t49 + t54) * qJD(2) + (-t41 - t92) * t91, t69, t21 * t34 + t33 * t22 - t29 * t32 - t31 * t64, t25, -t77, 0, qJD(2) * t93 + t46 * t22 + t37 * t32 - t82, -t83 - t46 * t21 + t37 * t31 + (-qJD(1) * t33 + t29) * qJD(2), t14 * t22 + t3 * t34 + t9 * t32 + t4 * t64 - t82, -t16 * t21 - t17 * t22 + t7 * t29 - t6 * t64 - t62, t14 * t21 - t4 * t29 + t3 * t33 - t9 * t31 + t83, t1 * t17 + t3 * t14 + t2 * t16 + t9 * t4 + t5 * t7 + t8 * t6; 0, 0, 0, 0, -t60, -t60 * qJ(2), -t60 * t55, -t60 * t56, 0, (-t49 - t91) * qJD(1), 0, 0, 0, 0, 0, t94, -t68, t94, -t34 * t22 - t32 * t64 - t69, t68, -t9 * qJD(1) + t62; 0, 0, 0, 0, 0, 0, 0, 0, -t84 * t60, t41 * t70 + t53, 0, 0, 0, 0, 0, t63, -t12, t63, -t26 - t90, t12, t8 * t64 + (-qJD(5) - t5) * t29 + t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, -t26 + t90, 0, t40 + (t29 - t74) * qJD(4), 0, -t37 * t29 + t66, t72 * qJD(4) - t20 + (qJD(3) + t37) * t64, -t15 * t64 + t66 - t89, pkin(4) * t21 - t22 * qJ(5) + (-t11 + t8) * t29 + (t5 - t76) * t64, t15 * t29 - t9 * t64 + (0.2e1 * qJD(5) - t72) * qJD(4) + t61, -t2 * pkin(4) + t1 * qJ(5) - t5 * t11 - t9 * t15 + t76 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, 0, -qJD(4) ^ 2 - t90, -t8 * qJD(4) + t2 + t89;];
tauc_reg = t13;
