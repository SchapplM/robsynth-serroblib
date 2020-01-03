% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% tauc_reg [5x21]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPRR7_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR7_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:59:53
% EndTime: 2019-12-31 17:59:55
% DurationCPUTime: 0.55s
% Computational Cost: add. (487->127), mult. (1062->200), div. (0->0), fcn. (619->6), ass. (0->79)
t36 = sin(qJ(4));
t61 = t36 * qJD(1);
t27 = qJD(5) + t61;
t92 = -qJD(5) + t27;
t38 = cos(qJ(4));
t35 = sin(qJ(5));
t65 = qJD(5) * t35;
t57 = t38 * t65;
t37 = cos(qJ(5));
t60 = t37 * qJD(4);
t32 = t38 ^ 2;
t69 = qJD(1) * t32;
t91 = -(-t27 * t36 + t69) * t60 + t27 * t57;
t68 = qJD(1) * t38;
t54 = t37 * t68;
t62 = t35 * qJD(4);
t22 = t54 + t62;
t55 = t36 * t62;
t9 = -qJD(1) * t55 + qJD(5) * t22;
t90 = 2 * qJD(3);
t26 = -cos(pkin(8)) * pkin(1) - pkin(2) - pkin(6);
t18 = t26 * qJD(1) + qJD(3);
t44 = t36 * qJD(2) - t38 * t18;
t4 = -qJD(4) * pkin(4) + t44;
t89 = t35 * t4;
t88 = t36 * t9;
t87 = t37 * t4;
t11 = t38 * qJD(2) + t36 * t18;
t7 = t11 * qJD(4);
t86 = t7 * t35;
t85 = t7 * t37;
t41 = -t36 * t60 - t57;
t8 = t41 * qJD(1) + qJD(5) * t60;
t84 = t8 * t35;
t66 = qJD(4) * t38;
t83 = t22 * t66 + t8 * t36;
t20 = t35 * t68 - t60;
t82 = t20 * t27;
t81 = t20 * t38;
t80 = t22 * t27;
t79 = t27 * t35;
t28 = sin(pkin(8)) * pkin(1) + qJ(3);
t16 = t36 * pkin(4) - t38 * pkin(7) + t28;
t13 = t16 * qJD(1);
t78 = t35 * t13;
t77 = t35 * t36;
t76 = t37 * t13;
t75 = t37 * t27;
t74 = t37 * t38;
t39 = qJD(4) ^ 2;
t73 = t39 * t36;
t72 = t39 * t38;
t71 = t36 ^ 2 - t32;
t40 = qJD(1) ^ 2;
t70 = -t39 - t40;
t24 = qJD(1) * t28;
t67 = qJD(4) * t36;
t64 = qJD(5) * t37;
t63 = t24 * qJD(1);
t59 = qJD(1) * qJD(4);
t58 = t35 * t69;
t56 = t27 * t64;
t47 = pkin(4) * t38 + pkin(7) * t36;
t19 = t47 * qJD(4) + qJD(3);
t14 = t19 * qJD(1);
t6 = t44 * qJD(4);
t53 = t37 * t14 + t35 * t6;
t5 = qJD(4) * pkin(7) + t11;
t52 = t26 * t27 + t5;
t51 = t38 * t59;
t50 = t27 + t61;
t48 = qJD(5) * t36 + qJD(1);
t1 = -t35 * t5 + t76;
t2 = t37 * t5 + t78;
t46 = -t35 * t14 + t37 * t6;
t43 = t27 * t55 - t38 * t56;
t42 = -pkin(7) * t66 + t36 * t4;
t23 = t47 * qJD(1);
t3 = [0, 0, 0, 0, 0, qJD(1) * t90, t24 * t90, -0.2e1 * t36 * t51, 0.2e1 * t71 * t59, -t73, -t72, 0, t24 * t66 - t26 * t73 + (t28 * t66 + t36 * t90) * qJD(1), -t24 * t67 - t26 * t72 + (-t28 * t67 + t38 * t90) * qJD(1), t41 * t22 + t8 * t74, (t20 * t37 + t22 * t35) * t67 + (-t84 - t37 * t9 + (t20 * t35 - t22 * t37) * qJD(5)) * t38, t83 - t91, -t88 + (-t58 - t81) * qJD(4) + t43, t50 * t66, (-t16 * t65 + t37 * t19) * t27 + ((t20 * t26 - t89) * qJD(4) + (-t52 * t37 - t78) * qJD(5) + t53) * t36 + (t4 * t64 - t26 * t9 + t86 + (-t26 * t79 + (t37 * t16 - t26 * t77) * qJD(1) + t1) * qJD(4)) * t38, -(t16 * t64 + t35 * t19) * t27 + ((t22 * t26 - t87) * qJD(4) + (t52 * t35 - t76) * qJD(5) + t46) * t36 + (-t4 * t65 - t26 * t8 + t85 + (-t26 * t75 - (t37 * t36 * t26 + t35 * t16) * qJD(1) - t2) * qJD(4)) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, t73, 0, 0, 0, 0, 0, t88 + (-t58 + t81) * qJD(4) + t43, t83 + t91; 0, 0, 0, 0, 0, -t40, -t63, 0, 0, 0, 0, 0, t70 * t36, t70 * t38, 0, 0, 0, 0, 0, -t38 * t9 - t48 * t75 + (-t50 * t38 * t35 + t20 * t36) * qJD(4), -t38 * t8 + t48 * t79 + (-t27 * t74 + (t22 - t54) * t36) * qJD(4); 0, 0, 0, 0, 0, 0, 0, t38 * t40 * t36, -t71 * t40, 0, 0, 0, -t38 * t63, t24 * t61, t22 * t75 + t84, (t8 - t82) * t37 + (-t9 - t80) * t35, t56 + (t36 * t75 + (-t22 + t62) * t38) * qJD(1), -t27 * t65 + (-t27 * t77 + (t20 + t60) * t38) * qJD(1), -t27 * t68, -pkin(4) * t9 - t85 - (t37 * t23 + t35 * t44) * t27 - t11 * t20 + (-pkin(7) * t75 + t89) * qJD(5) + (-t1 * t38 + t42 * t35) * qJD(1), -pkin(4) * t8 + t86 + (t35 * t23 - t37 * t44) * t27 - t11 * t22 + (pkin(7) * t79 + t87) * qJD(5) + (t2 * t38 + t42 * t37) * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22 * t20, -t20 ^ 2 + t22 ^ 2, t8 + t82, t80 - t9, t51, t92 * t2 - t4 * t22 + t53, t92 * t1 + t4 * t20 + t46;];
tauc_reg = t3;
