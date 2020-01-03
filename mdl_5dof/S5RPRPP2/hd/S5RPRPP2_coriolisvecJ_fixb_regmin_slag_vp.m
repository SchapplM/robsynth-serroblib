% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
% 
% Output:
% tauc_reg [5x19]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPP2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:11:12
% EndTime: 2019-12-31 18:11:13
% DurationCPUTime: 0.46s
% Computational Cost: add. (487->117), mult. (1095->159), div. (0->0), fcn. (517->4), ass. (0->84)
t93 = pkin(3) + pkin(4);
t70 = t93 * qJD(3);
t35 = sin(pkin(7)) * pkin(1) + pkin(6);
t30 = t35 * qJD(1);
t50 = sin(qJ(3));
t27 = t50 * t30;
t51 = cos(qJ(3));
t19 = t51 * qJD(2) - t27;
t76 = qJ(5) * qJD(1);
t8 = t50 * t76 + t19;
t85 = qJD(4) - t8;
t3 = -t70 + t85;
t96 = qJD(4) - t19;
t45 = qJD(3) * qJ(4);
t20 = t50 * qJD(2) + t51 * t30;
t9 = -t51 * t76 + t20;
t6 = t45 + t9;
t36 = -cos(pkin(7)) * pkin(1) - pkin(2);
t58 = t50 * qJ(4) - t36;
t18 = t93 * t51 + t58;
t83 = qJD(1) * t18;
t4 = qJD(5) + t83;
t84 = qJD(5) + t4;
t95 = t84 * t50;
t14 = t45 + t20;
t10 = -qJD(3) * pkin(3) + t96;
t94 = (t19 + t27) * qJD(3);
t92 = t51 * t6;
t53 = qJD(3) ^ 2;
t91 = t53 * t50;
t42 = t53 * t51;
t73 = qJD(2) * qJD(3);
t80 = qJD(3) * t51;
t17 = t30 * t80 + t50 * t73;
t74 = qJD(1) * qJD(3);
t68 = t51 * t74;
t79 = t50 * qJD(4);
t90 = qJ(4) * t68 + qJD(1) * t79;
t39 = t51 * t73;
t44 = qJD(3) * qJD(4);
t89 = t39 + 0.2e1 * t44;
t46 = t50 ^ 2;
t47 = t51 ^ 2;
t88 = t46 - t47;
t87 = qJ(4) * t51;
t86 = qJ(5) - t35;
t23 = -t51 * pkin(3) - t58;
t15 = qJD(1) * t23;
t31 = qJD(1) * t36;
t82 = qJD(1) * t50;
t81 = qJD(3) * t50;
t78 = t50 * qJD(5);
t77 = t51 * qJD(5);
t75 = qJ(5) * qJD(3);
t54 = qJD(1) ^ 2;
t72 = t50 * t54 * t51;
t71 = t51 * t75;
t69 = t50 * t74;
t25 = t86 * t51;
t56 = -t93 * t50 + t87;
t11 = t56 * qJD(3) + t79;
t5 = -t93 * t69 + t90;
t66 = qJD(1) * t11 + t5;
t65 = t4 + t83;
t64 = 0.2e1 * t15;
t62 = 0.2e1 * t68;
t61 = t20 * qJD(3) - t17;
t60 = pkin(3) * t50 - t87;
t59 = 0.2e1 * qJD(3) * t31;
t7 = -t30 * t81 + t39 + t44;
t16 = pkin(3) * t69 - t90;
t22 = t60 * qJD(3) - t79;
t57 = -qJD(1) * t22 - t35 * t53 - t16;
t55 = t17 * t50 + t7 * t51 + (t10 * t51 - t14 * t50) * qJD(3);
t33 = qJ(5) * t69;
t32 = -t46 * t54 - t53;
t29 = t60 * qJD(1);
t24 = t86 * t50;
t21 = t56 * qJD(1);
t13 = -qJD(3) * t25 - t78;
t12 = t86 * t81 - t77;
t2 = (-t71 - t78) * qJD(1) + t17;
t1 = -qJD(1) * t77 + t33 + t7;
t26 = [0, 0, 0, 0, t50 * t62, -0.2e1 * t88 * t74, t42, -t91, 0, -t35 * t42 + t50 * t59, t35 * t91 + t51 * t59, t57 * t51 + t64 * t81, t55, t57 * t50 - t64 * t80, t15 * t22 + t16 * t23 + t55 * t35, t66 * t51 + (-t65 * t50 - t13) * qJD(3), t66 * t50 + (t65 * t51 + t12) * qJD(3), -t1 * t51 - t2 * t50 + (-t3 * t51 + t50 * t6) * qJD(3) + (-t12 * t51 - t13 * t50 + (t24 * t51 - t25 * t50) * qJD(3)) * qJD(1), -t1 * t25 + t4 * t11 + t6 * t12 + t3 * t13 + t5 * t18 - t2 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91, -t42, -t91, 0, t42, -t17 * t51 + t7 * t50 + (t10 * t50 + t14 * t51) * qJD(3), -t91, t42, 0, t1 * t50 - t2 * t51 + (t3 * t50 + t92) * qJD(3); 0, 0, 0, 0, -t72, t88 * t54, 0, 0, 0, -t31 * t82 + t61, -t31 * t51 * qJD(1) - t39 + t94, (-t15 * t50 + t29 * t51) * qJD(1) + t61, 0, -t94 + (t15 * t51 + t29 * t50) * qJD(1) + t89, -t17 * pkin(3) + t7 * qJ(4) - t10 * t20 + t96 * t14 - t15 * t29, t9 * qJD(3) + ((-t21 + t75) * t51 + t95) * qJD(1) - t17, t33 + (-t8 - t27) * qJD(3) + (-t21 * t50 - t84 * t51) * qJD(1) + t89, 0, t1 * qJ(4) - t2 * t93 - t4 * t21 - t3 * t9 + t85 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, 0, t32, -t14 * qJD(3) + t15 * t82 + t17, -t72, t32, 0, -t6 * qJD(3) + (-t71 - t95) * qJD(1) + t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t69, t62, (-t46 - t47) * t54, (t92 + (t3 - t70) * t50) * qJD(1) + t90;];
tauc_reg = t26;
