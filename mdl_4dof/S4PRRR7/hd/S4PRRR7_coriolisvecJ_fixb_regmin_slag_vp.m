% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4PRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% tauc_reg [4x18]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PRRR7_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR7_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:36:48
% EndTime: 2019-12-31 16:36:50
% DurationCPUTime: 0.70s
% Computational Cost: add. (401->134), mult. (1140->231), div. (0->0), fcn. (808->8), ass. (0->88)
t38 = cos(qJ(3));
t69 = t38 * qJD(2);
t27 = -qJD(4) + t69;
t96 = -qJD(4) - t27;
t34 = sin(qJ(4));
t35 = sin(qJ(3));
t37 = cos(qJ(4));
t72 = qJD(4) * t37;
t61 = t35 * t72;
t67 = qJD(3) * qJD(4);
t71 = t34 * qJD(3);
t11 = qJD(2) * (t38 * t71 + t61) + t34 * t67;
t36 = sin(qJ(2));
t32 = sin(pkin(4));
t79 = qJD(1) * t32;
t60 = t36 * t79;
t23 = qJD(2) * pkin(6) + t60;
t33 = cos(pkin(4));
t39 = cos(qJ(2));
t77 = qJD(2) * t32;
t64 = t39 * t77;
t44 = qJD(1) * (qJD(3) * t33 + t64);
t74 = qJD(3) * t38;
t4 = t23 * t74 + t35 * t44;
t95 = t4 * t34;
t94 = t4 * t37;
t70 = t37 * qJD(3);
t57 = t38 * t70;
t76 = qJD(2) * t35;
t58 = t34 * t76;
t10 = qJD(2) * t57 - qJD(4) * t58 + t37 * t67;
t93 = t10 * t34;
t18 = t58 - t70;
t92 = t18 * t27;
t20 = t37 * t76 + t71;
t91 = t20 * t27;
t90 = t27 * t34;
t89 = t32 * t36;
t88 = t32 * t39;
t41 = qJD(2) ^ 2;
t87 = t32 * t41;
t86 = t34 * t38;
t85 = t37 * t27;
t84 = t37 * t38;
t40 = qJD(3) ^ 2;
t83 = t40 * t35;
t82 = t40 * t38;
t30 = t35 ^ 2;
t81 = -t38 ^ 2 + t30;
t80 = qJD(2) * pkin(2);
t78 = qJD(1) * t33;
t75 = qJD(3) * t35;
t73 = qJD(4) * t34;
t68 = qJD(2) * qJD(3);
t66 = t36 * t87;
t65 = t36 * t77;
t63 = t35 * t73;
t62 = t27 * t72;
t59 = t39 * t79;
t55 = t35 * t68;
t54 = t35 * t64;
t53 = t38 * t64;
t52 = pkin(3) * t35 - pkin(7) * t38;
t22 = t52 * qJD(3);
t51 = -t22 + t60;
t25 = -t38 * pkin(3) - t35 * pkin(7) - pkin(2);
t15 = t25 * qJD(2) - t59;
t14 = t38 * t23 + t35 * t78;
t9 = qJD(3) * pkin(7) + t14;
t1 = t37 * t15 - t34 * t9;
t2 = t34 * t15 + t37 * t9;
t49 = qJD(2) * t30 - t27 * t38;
t17 = t33 * t35 + t38 * t89;
t48 = -t17 * t34 - t37 * t88;
t47 = -t17 * t37 + t34 * t88;
t16 = -t33 * t38 + t35 * t89;
t13 = -t35 * t23 + t38 * t78;
t46 = qJD(2) * t80;
t43 = -0.2e1 * qJD(3) * t80;
t3 = -t23 * t75 + t38 * t44;
t8 = -qJD(3) * pkin(3) - t13;
t42 = qJD(3) * t8 + qJD(4) * t15 - t27 * t59 + t3;
t21 = t52 * qJD(2);
t12 = (t22 + t60) * qJD(2);
t7 = t37 * t12;
t6 = t17 * qJD(3) + t54;
t5 = -qJD(3) * t16 + t53;
t19 = [0, 0, -t66, -t39 * t87, 0, 0, 0, 0, 0, -t38 * t66 + (-t6 - t54) * qJD(3), t35 * t66 + (-t5 - t53) * qJD(3), 0, 0, 0, 0, 0, -(t47 * qJD(4) - t5 * t34 + t37 * t65) * t27 + t48 * t55 + t6 * t18 + t16 * t11, (t48 * qJD(4) + t34 * t65 + t5 * t37) * t27 + t47 * t55 + t6 * t20 + t16 * t10; 0, 0, 0, 0, 0.2e1 * t38 * t55, -0.2e1 * t81 * t68, t82, -t83, 0, -pkin(6) * t82 + t35 * t43, pkin(6) * t83 + t38 * t43, t10 * t37 * t35 + (t57 - t63) * t20, (-t18 * t37 - t20 * t34) * t74 + (-t93 - t11 * t37 + (t18 * t34 - t20 * t37) * qJD(4)) * t35, t27 * t63 - t10 * t38 + (t20 * t35 + t49 * t37) * qJD(3), t27 * t61 + t11 * t38 + (-t18 * t35 - t49 * t34) * qJD(3), (-t27 - t69) * t75, (t25 * t73 + t51 * t37) * t27 + (t9 * t72 - t7 + (qJD(3) * t18 + t62) * pkin(6) + t42 * t34) * t38 + (-t18 * t59 + t8 * t72 + pkin(6) * t11 + t95 + (-pkin(6) * t90 + (-pkin(6) * t86 + t37 * t25) * qJD(2) + t1) * qJD(3)) * t35, (t25 * t72 - t51 * t34) * t27 + (qJD(3) * pkin(6) * t20 + (t12 + (-pkin(6) * t27 - t9) * qJD(4)) * t34 + t42 * t37) * t38 + (-t20 * t59 - t8 * t73 + pkin(6) * t10 + t94 + (-pkin(6) * t85 - (pkin(6) * t84 + t34 * t25) * qJD(2) - t2) * qJD(3)) * t35; 0, 0, 0, 0, -t35 * t41 * t38, t81 * t41, 0, 0, 0, t35 * t46, t38 * t46, -t20 * t85 + t93, (t10 + t92) * t37 + (-t11 + t91) * t34, -t62 + (t27 * t84 + (-t20 + t71) * t35) * qJD(2), t27 * t73 + (-t27 * t86 + (t18 + t70) * t35) * qJD(2), t27 * t76, -pkin(3) * t11 - t94 + (-t34 * t13 + t37 * t21) * t27 - t14 * t18 + (pkin(7) * t85 + t8 * t34) * qJD(4) + (-t1 * t35 + (-pkin(7) * t75 - t38 * t8) * t34) * qJD(2), -pkin(3) * t10 + t95 - (t37 * t13 + t34 * t21) * t27 - t14 * t20 + (-pkin(7) * t90 + t8 * t37) * qJD(4) + (-t8 * t84 + (-pkin(7) * t70 + t2) * t35) * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20 * t18, -t18 ^ 2 + t20 ^ 2, t10 - t92, -t11 - t91, t55, t96 * t2 - t8 * t20 - t34 * t3 + t7, t96 * t1 - t34 * t12 + t8 * t18 - t37 * t3;];
tauc_reg = t19;
