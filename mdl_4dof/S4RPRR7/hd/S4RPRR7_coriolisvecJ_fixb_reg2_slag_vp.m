% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RPRR7
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
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPRR7_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR7_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR7_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:54:04
% EndTime: 2019-12-31 16:54:07
% DurationCPUTime: 0.76s
% Computational Cost: add. (1638->176), mult. (4449->252), div. (0->0), fcn. (3198->6), ass. (0->94)
t110 = cos(qJ(3));
t64 = sin(pkin(7));
t65 = cos(pkin(7));
t67 = sin(qJ(3));
t49 = t110 * t64 + t67 * t65;
t43 = t49 * qJD(1);
t66 = sin(qJ(4));
t68 = cos(qJ(4));
t34 = t66 * qJD(3) + t68 * t43;
t86 = t110 * t65;
t56 = qJD(1) * t86;
t96 = t67 * t64;
t87 = qJD(1) * t96;
t41 = -t56 + t87;
t39 = qJD(4) + t41;
t84 = t39 * t66;
t115 = t34 * t84;
t73 = t86 - t96;
t70 = t73 * qJD(2);
t95 = pkin(5) + qJ(2);
t53 = t95 * t64;
t50 = qJD(1) * t53;
t54 = t95 * t65;
t51 = qJD(1) * t54;
t75 = t110 * t50 + t67 * t51;
t11 = qJD(1) * t70 - t75 * qJD(3);
t55 = qJD(3) * t56;
t37 = qJD(3) * t87 - t55;
t46 = t49 * qJD(3);
t38 = qJD(1) * t46;
t20 = t38 * pkin(3) + t37 * pkin(6);
t59 = -t65 * pkin(2) - pkin(1);
t52 = t59 * qJD(1) + qJD(2);
t16 = t41 * pkin(3) - t43 * pkin(6) + t52;
t28 = t110 * t51 - t67 * t50;
t23 = qJD(3) * pkin(6) + t28;
t5 = t68 * t16 - t66 * t23;
t1 = qJD(4) * t5 + t68 * t11 + t66 * t20;
t81 = -t5 * t39 + t1;
t6 = t66 * t16 + t68 * t23;
t2 = -qJD(4) * t6 - t66 * t11 + t68 * t20;
t114 = t6 * t39 + t2;
t92 = qJD(4) * t34;
t15 = -t66 * t37 + t92;
t113 = t43 ^ 2;
t71 = t49 * qJD(2);
t12 = qJD(1) * t71 + t28 * qJD(3);
t74 = -t110 * t53 - t67 * t54;
t109 = t12 * t74;
t108 = t12 * t66;
t89 = t68 * qJD(3);
t91 = qJD(4) * t66;
t14 = -qJD(4) * t89 + t68 * t37 + t43 * t91;
t107 = t14 * t66;
t106 = t15 * t68;
t32 = t66 * t43 - t89;
t105 = t32 * t41;
t104 = t34 * t32;
t103 = t34 * t43;
t102 = t38 * t73;
t101 = t43 * t32;
t100 = t43 * t41;
t99 = t49 * t68;
t13 = t66 * t15;
t97 = t66 * t38;
t36 = t68 * t38;
t90 = qJD(4) * t68;
t94 = -t32 * t90 - t13;
t93 = t64 ^ 2 + t65 ^ 2;
t88 = qJD(1) * qJD(2);
t85 = t93 * qJD(1) ^ 2;
t83 = t39 * t68;
t82 = t5 * t68 + t6 * t66;
t26 = -pkin(3) * t73 - t49 * pkin(6) + t59;
t31 = t110 * t54 - t67 * t53;
t9 = t68 * t26 - t66 * t31;
t10 = t66 * t26 + t68 * t31;
t79 = 0.2e1 * t93 * t88;
t78 = -t39 * t91 - t41 * t84 + t36;
t45 = t73 * qJD(3);
t77 = t45 * t66 + t49 * t90;
t76 = t45 * t68 - t49 * t91;
t22 = -qJD(3) * pkin(3) + t75;
t72 = -pkin(6) * t38 + t39 * t22;
t40 = t41 ^ 2;
t25 = t46 * pkin(3) - t45 * pkin(6);
t24 = t43 * pkin(3) + t41 * pkin(6);
t18 = t31 * qJD(3) + t71;
t17 = t74 * qJD(3) + t70;
t8 = t66 * t24 - t68 * t75;
t7 = t68 * t24 + t66 * t75;
t4 = -t10 * qJD(4) - t66 * t17 + t68 * t25;
t3 = t9 * qJD(4) + t68 * t17 + t66 * t25;
t19 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, qJ(2) * t79, -t37 * t49 + t43 * t45, -t37 * t73 - t49 * t38 - t45 * t41 - t43 * t46, t45 * qJD(3), t41 * t46 - t102, -t46 * qJD(3), 0, -t18 * qJD(3) + t59 * t38 + t52 * t46, -t17 * qJD(3) - t59 * t37 + t52 * t45, t11 * t73 + t12 * t49 - t17 * t41 + t18 * t43 - t28 * t46 - t31 * t38 + t37 * t74 + t45 * t75, t11 * t31 + t28 * t17 + t18 * t75 - t109, -t14 * t99 + t76 * t34, -(t32 * t68 + t34 * t66) * t45 + (t107 - t106 + (t32 * t66 - t34 * t68) * qJD(4)) * t49, t14 * t73 + t34 * t46 + t49 * t36 + t76 * t39, t49 * t13 + t77 * t32, t15 * t73 - t32 * t46 - t77 * t39 - t49 * t97, t39 * t46 - t102, t49 * t108 - t15 * t74 + t18 * t32 - t2 * t73 + t77 * t22 + t9 * t38 + t4 * t39 + t5 * t46, t1 * t73 - t10 * t38 + t12 * t99 + t14 * t74 + t18 * t34 + t76 * t22 - t3 * t39 - t6 * t46, -t10 * t15 + t9 * t14 - t3 * t32 - t4 * t34 - t82 * t45 + (-t1 * t66 - t2 * t68 + (t5 * t66 - t6 * t68) * qJD(4)) * t49, t1 * t10 + t22 * t18 + t2 * t9 + t6 * t3 + t5 * t4 - t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85, -qJ(2) * t85, 0, 0, 0, 0, 0, 0, 0.2e1 * t43 * qJD(3), t55 + (-t41 - t87) * qJD(3), -t40 - t113, t28 * t41 - t43 * t75, 0, 0, 0, 0, 0, 0, t78 - t101, -t39 ^ 2 * t68 - t103 - t97, (t14 - t105) * t68 + t115 + t94, t114 * t68 - t22 * t43 + t81 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, -t40 + t113, t55 + (t41 - t87) * qJD(3), -t100, 0, 0, -(qJD(2) + t52) * t43, t52 * t41 - t73 * t88, 0, 0, t34 * t83 - t107, (-t14 - t105) * t68 - t115 + t94, t39 * t83 - t103 + t97, t32 * t84 - t106, t78 + t101, -t39 * t43, -pkin(3) * t15 - t12 * t68 - t28 * t32 - t5 * t43 + (-pkin(6) * t90 - t7) * t39 + t72 * t66, pkin(3) * t14 + t108 - t28 * t34 + t6 * t43 + (pkin(6) * t91 + t8) * t39 + t72 * t68, t8 * t32 + t7 * t34 + ((-t15 + t92) * pkin(6) + t81) * t68 + ((qJD(4) * t32 - t14) * pkin(6) - t114) * t66, -t12 * pkin(3) - t22 * t28 - t5 * t7 - t6 * t8 + (-t82 * qJD(4) + t1 * t68 - t2 * t66) * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104, -t32 ^ 2 + t34 ^ 2, t32 * t39 - t14, -t104, t34 * t39 - t15, t38, -t22 * t34 + t114, t22 * t32 - t81, 0, 0;];
tauc_reg = t19;
