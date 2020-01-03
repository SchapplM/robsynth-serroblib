% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPRR8_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR8_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR8_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR8_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:55:15
% EndTime: 2019-12-31 16:55:17
% DurationCPUTime: 0.41s
% Computational Cost: add. (774->110), mult. (1631->164), div. (0->0), fcn. (953->4), ass. (0->84)
t74 = 2 * qJD(1);
t55 = (-pkin(1) - pkin(5));
t40 = t55 * qJD(1) + qJD(2);
t54 = cos(qJ(3));
t82 = qJD(1) * t54;
t25 = -pkin(6) * t82 + t54 * t40;
t21 = qJD(3) * pkin(3) + t25;
t65 = pkin(6) * qJD(1) - t40;
t80 = qJD(3) * t54;
t23 = t65 * t80;
t53 = cos(qJ(4));
t97 = (qJD(4) * t21 - t23) * t53;
t51 = sin(qJ(4));
t52 = sin(qJ(3));
t31 = t51 * t54 + t53 * t52;
t26 = t31 * qJD(1);
t48 = qJD(3) + qJD(4);
t96 = pkin(6) - t55;
t90 = t53 * t54;
t62 = t48 * t90;
t79 = qJD(4) * t51;
t81 = qJD(3) * t52;
t16 = -t51 * t81 - t52 * t79 + t62;
t95 = t16 * t48;
t71 = t53 * t82;
t83 = qJD(1) * t52;
t72 = t51 * t83;
t28 = t71 - t72;
t94 = t28 * t26;
t45 = t52 * pkin(3) + qJ(2);
t37 = t45 * qJD(1);
t93 = t37 * t28;
t24 = -pkin(6) * t83 + t52 * t40;
t92 = t51 * t24;
t91 = t53 * t24;
t56 = qJD(3) ^ 2;
t89 = t56 * t52;
t88 = t56 * t54;
t76 = qJD(1) * qJD(3);
t68 = t52 * t76;
t87 = -qJD(4) * t72 - t51 * t68;
t86 = t52 ^ 2 - t54 ^ 2;
t57 = qJD(1) ^ 2;
t85 = -t56 - t57;
t84 = t57 * qJ(2);
t78 = t37 * qJD(1);
t77 = qJ(2) * qJD(3);
t75 = t54 * t57 * t52;
t73 = pkin(3) * t82;
t70 = qJD(2) * t74;
t36 = t96 * t54;
t69 = -pkin(3) * t48 - t21;
t22 = t65 * t81;
t67 = t53 * t22 + t51 * t23;
t66 = t51 * t22 - t24 * t79;
t63 = t54 * t68;
t41 = pkin(3) * t80 + qJD(2);
t15 = t48 * t31;
t12 = t15 * qJD(1);
t32 = -t51 * t52 + t90;
t61 = -t32 * t12 - t28 * t15;
t13 = qJD(1) * t62 + t87;
t60 = t31 * t13 + t16 * t26;
t9 = t51 * t21 + t91;
t35 = t96 * t52;
t18 = -t53 * t35 - t51 * t36;
t17 = t51 * t35 - t53 * t36;
t59 = t37 * t26 - t66;
t1 = t66 + t97;
t2 = -qJD(4) * t9 + t67;
t8 = t53 * t21 - t92;
t58 = t1 * t31 - t8 * t15 + t9 * t16 + t2 * t32;
t47 = qJ(2) * t70;
t34 = t41 * qJD(1);
t30 = qJD(3) * t36;
t29 = t96 * t81;
t14 = t15 * t48;
t11 = t53 * t25 - t92;
t10 = -t51 * t25 - t91;
t7 = -t26 ^ 2 + t28 ^ 2;
t6 = -t87 + (t28 - t71) * t48;
t4 = -t18 * qJD(4) + t53 * t29 + t51 * t30;
t3 = t17 * qJD(4) + t51 * t29 - t53 * t30;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, t47, -0.2e1 * t63, 0.2e1 * t86 * t76, -t89, 0.2e1 * t63, -t88, 0, -t55 * t89 + (qJD(2) * t52 + t54 * t77) * t74, -t55 * t88 + (qJD(2) * t54 - t52 * t77) * t74, 0, t47, t61, t12 * t31 - t32 * t13 + t15 * t26 - t28 * t16, -t14, t60, -t95, 0, t45 * t13 + t37 * t16 + t41 * t26 + t34 * t31 + t4 * t48, -t45 * t12 - t37 * t15 + t41 * t28 - t3 * t48 + t34 * t32, t17 * t12 - t18 * t13 - t3 * t26 - t4 * t28 - t58, t1 * t18 + t2 * t17 + t9 * t3 + t34 * t45 + t37 * t41 + t8 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, -t84, 0, 0, 0, 0, 0, 0, t85 * t52, t85 * t54, 0, -t84, 0, 0, 0, 0, 0, 0, -qJD(1) * t26 - t14, -qJD(1) * t28 - t95, -t60 - t61, t58 - t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, -t86 * t57, 0, -t75, 0, 0, -t54 * t84, t52 * t84, 0, 0, t94, t7, 0, -t94, t6, 0, -t26 * t73 - t10 * t48 - t93 + (t69 * t51 - t91) * qJD(4) + t67, -t28 * t73 + t11 * t48 + (t69 * qJD(4) + t23) * t53 + t59, (t10 + t9) * t28 + (t11 - t8) * t26 + (t12 * t53 - t13 * t51 + (-t26 * t53 + t28 * t51) * qJD(4)) * pkin(3), -t8 * t10 - t9 * t11 + (-t54 * t78 + t1 * t51 + t2 * t53 + (-t51 * t8 + t53 * t9) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, t7, 0, -t94, t6, 0, t9 * t48 + t2 - t93, t8 * t48 + t59 - t97, 0, 0;];
tauc_reg = t5;
