% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RPRR4
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
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPRR4_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR4_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:50:35
% EndTime: 2019-12-31 16:50:38
% DurationCPUTime: 0.74s
% Computational Cost: add. (900->149), mult. (2224->233), div. (0->0), fcn. (1299->6), ass. (0->94)
t35 = sin(pkin(7)) * pkin(1) + pkin(5);
t31 = t35 * qJD(1);
t44 = sin(qJ(3));
t46 = cos(qJ(3));
t20 = t46 * qJD(2) - t44 * t31;
t16 = t20 * qJD(3);
t59 = pkin(3) * t44 - pkin(6) * t46;
t30 = t59 * qJD(3);
t24 = qJD(1) * t30;
t43 = sin(qJ(4));
t45 = cos(qJ(4));
t21 = t44 * qJD(2) + t46 * t31;
t15 = qJD(3) * pkin(6) + t21;
t36 = -cos(pkin(7)) * pkin(1) - pkin(2);
t25 = -t46 * pkin(3) - t44 * pkin(6) + t36;
t18 = t25 * qJD(1);
t56 = t43 * t15 - t45 * t18;
t1 = -t56 * qJD(4) + t45 * t16 + t43 * t24;
t78 = t46 * qJD(1);
t34 = -qJD(4) + t78;
t110 = -t56 * t34 + t1;
t6 = t45 * t15 + t43 * t18;
t2 = -qJD(4) * t6 - t43 * t16 + t45 * t24;
t109 = t6 * t34 - t2;
t39 = t44 ^ 2;
t54 = qJD(1) * t39 - t34 * t46;
t82 = qJD(4) * t43;
t71 = t44 * t82;
t79 = t45 * qJD(3);
t108 = -t34 * t71 - t54 * t79;
t81 = qJD(4) * t45;
t70 = t44 * t81;
t80 = t43 * qJD(3);
t51 = t46 * t80 + t70;
t76 = qJD(3) * qJD(4);
t13 = t51 * qJD(1) + t43 * t76;
t14 = -qJD(3) * pkin(3) - t20;
t105 = t14 * t43;
t104 = t14 * t45;
t17 = t21 * qJD(3);
t103 = t17 * t43;
t102 = t17 * t44;
t101 = t17 * t45;
t86 = qJD(1) * t44;
t68 = t43 * t86;
t26 = t68 - t79;
t100 = t26 * t34;
t99 = t26 * t44;
t28 = t45 * t86 + t80;
t98 = t28 * t26;
t97 = t28 * t34;
t96 = t34 * t43;
t95 = t34 * t45;
t94 = t43 * t46;
t93 = t44 * t45;
t92 = t45 * t46;
t91 = t46 * t13;
t47 = qJD(3) ^ 2;
t90 = t47 * t44;
t89 = t47 * t46;
t72 = t46 * t79;
t88 = -t13 * t93 - t26 * t72;
t87 = -t46 ^ 2 + t39;
t32 = qJD(1) * t36;
t85 = qJD(3) * t44;
t84 = qJD(3) * t46;
t83 = qJD(4) * t26;
t77 = qJD(1) * qJD(3);
t48 = qJD(1) ^ 2;
t75 = t44 * t48 * t46;
t74 = t28 * t84;
t73 = t35 * t85;
t69 = t34 * t81;
t66 = t44 * t77;
t12 = -qJD(1) * t72 + qJD(4) * t68 - t45 * t76;
t65 = t12 * t46 + t28 * t85;
t64 = -t12 + t83;
t62 = t44 * t69;
t61 = t28 * t70;
t60 = t46 * t66;
t58 = -t43 * t6 + t45 * t56;
t57 = -t43 * t56 - t45 * t6;
t53 = 0.2e1 * qJD(3) * t32;
t10 = t43 * t25 + t35 * t92;
t9 = t45 * t25 - t35 * t94;
t52 = t54 * t43;
t50 = t58 * qJD(4) + t1 * t45 - t2 * t43;
t49 = t16 * t46 + t102 + (-t20 * t46 - t21 * t44) * qJD(3);
t29 = t59 * qJD(1);
t8 = t45 * t20 + t43 * t29;
t7 = -t43 * t20 + t45 * t29;
t4 = -t10 * qJD(4) + t45 * t30 + t43 * t73;
t3 = t9 * qJD(4) + t43 * t30 - t45 * t73;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t60, -0.2e1 * t87 * t77, t89, -0.2e1 * t60, -t90, 0, -t35 * t89 + t44 * t53, t35 * t90 + t46 * t53, t49, t49 * t35, -t12 * t93 + (-t71 + t72) * t28, -t61 + (-t74 + (t12 + t83) * t44) * t43 + t88, -t108 + t65, t13 * t43 * t44 + t51 * t26, t62 + t91 + (-t52 - t99) * qJD(3), (-t34 - t78) * t85, -t4 * t34 + (-t2 + (t26 * t35 + t105) * qJD(3)) * t46 + (t14 * t81 + t13 * t35 + t103 + (qJD(1) * t9 - t56) * qJD(3)) * t44, t3 * t34 + (t1 + (t28 * t35 + t104) * qJD(3)) * t46 + (-t14 * t82 - t12 * t35 + t101 + (-qJD(1) * t10 - t6) * qJD(3)) * t44, -t10 * t13 + t9 * t12 - t3 * t26 - t4 * t28 + t58 * t84 + (t57 * qJD(4) - t1 * t43 - t2 * t45) * t44, t1 * t10 + t2 * t9 + t6 * t3 - t56 * t4 + (t14 * t84 + t102) * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t90, -t89, 0, t16 * t44 - t17 * t46 + (-t20 * t44 + t21 * t46) * qJD(3), 0, 0, 0, 0, 0, 0, t62 - t91 + (-t52 + t99) * qJD(3), t108 + t65, t61 + (t64 * t44 + t74) * t43 + t88, (-t57 * qJD(3) - t17) * t46 + (qJD(3) * t14 + t50) * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75, t87 * t48, 0, t75, 0, 0, -t32 * t86, -t32 * t78, 0, 0, -t12 * t43 - t28 * t95, (-t12 + t100) * t45 + (-t13 + t97) * t43, -t69 + (t34 * t92 + (-t28 + t80) * t44) * qJD(1), -t13 * t45 - t26 * t96, t34 * t82 + (-t34 * t94 + (t26 + t79) * t44) * qJD(1), t34 * t86, -pkin(3) * t13 - t101 - t21 * t26 + t7 * t34 + (pkin(6) * t95 + t105) * qJD(4) + (t44 * t56 + (-pkin(6) * t85 - t14 * t46) * t43) * qJD(1), pkin(3) * t12 + t103 - t21 * t28 - t8 * t34 + (-pkin(6) * t96 + t104) * qJD(4) + (-t14 * t92 + (-pkin(6) * t79 + t6) * t44) * qJD(1), t8 * t26 + t7 * t28 + ((qJD(4) * t28 - t13) * pkin(6) + t110) * t45 + (t64 * pkin(6) + t109) * t43, -t17 * pkin(3) + pkin(6) * t50 - t14 * t21 + t56 * t7 - t6 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, -t26 ^ 2 + t28 ^ 2, -t12 - t100, -t98, -t13 - t97, t66, -t14 * t28 - t109, t14 * t26 - t110, 0, 0;];
tauc_reg = t5;
