% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% tauc_reg [5x19]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRPRP1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:28:40
% EndTime: 2019-12-05 15:28:42
% DurationCPUTime: 0.34s
% Computational Cost: add. (562->111), mult. (1447->145), div. (0->0), fcn. (1018->4), ass. (0->70)
t56 = sin(pkin(8));
t57 = cos(pkin(8));
t58 = sin(qJ(4));
t59 = cos(qJ(4));
t37 = t59 * t56 + t58 * t57;
t32 = t37 * qJD(2);
t89 = t32 ^ 2;
t79 = qJD(2) * t57;
t72 = t59 * t79;
t85 = t58 * t56;
t73 = qJD(2) * t85;
t30 = -t72 + t73;
t51 = -t57 * pkin(3) - pkin(2);
t40 = t51 * qJD(2) + qJD(3);
t6 = t30 * pkin(4) - t32 * qJ(5) + t40;
t88 = t6 * t32;
t87 = t32 * t30;
t75 = qJ(3) * qJD(2);
t39 = t56 * qJD(1) + t57 * t75;
t28 = pkin(6) * t79 + t39;
t86 = t58 * t28;
t84 = pkin(6) + qJ(3);
t35 = t37 * qJD(4);
t24 = qJD(2) * t35;
t77 = qJD(4) * t59;
t78 = qJD(4) * t58;
t34 = t56 * t78 - t57 * t77;
t83 = -t37 * t24 + t34 * t30;
t82 = t56 ^ 2 + t57 ^ 2;
t36 = -t59 * t57 + t85;
t41 = t84 * t56;
t42 = t84 * t57;
t64 = -t59 * t41 - t58 * t42;
t8 = -t36 * qJD(3) + t64 * qJD(4);
t81 = t8 * qJD(4);
t19 = -t58 * t41 + t59 * t42;
t9 = t37 * qJD(3) + t19 * qJD(4);
t80 = t9 * qJD(4);
t26 = t35 * qJD(4);
t53 = t57 * qJD(1);
t27 = -qJD(2) * t41 + t53;
t10 = t59 * t27 - t86;
t76 = qJD(5) - t10;
t74 = qJD(2) * qJD(3);
t71 = t58 * t74;
t70 = t59 * t74;
t69 = t10 + t86;
t68 = t82 * qJD(2);
t2 = t27 * t78 + t28 * t77 + t56 * t70 + t57 * t71;
t46 = qJD(4) * t72;
t23 = qJD(4) * t73 - t46;
t67 = t24 * pkin(4) + t23 * qJ(5);
t66 = -t36 * t23 + t35 * t32;
t11 = t58 * t27 + t59 * t28;
t65 = (-t56 * t75 + t53) * t56 - t39 * t57;
t63 = t11 * qJD(4) - t2;
t62 = -t27 * t77 + t56 * t71 - t57 * t70;
t61 = 0.2e1 * t32 * qJD(4);
t29 = t30 ^ 2;
t25 = t34 * qJD(4);
t15 = t36 * pkin(4) - t37 * qJ(5) + t51;
t14 = t32 * pkin(4) + t30 * qJ(5);
t13 = t46 + (t30 - t73) * qJD(4);
t12 = -t46 + (t30 + t73) * qJD(4);
t7 = qJD(4) * qJ(5) + t11;
t5 = -qJD(4) * pkin(4) + t76;
t4 = t35 * pkin(4) + t34 * qJ(5) - t37 * qJD(5);
t3 = -t32 * qJD(5) + t67;
t1 = (qJD(5) - t86) * qJD(4) - t62;
t16 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, t25, -t26, t66 + t83, -t25, t1 * t37 + t2 * t36 - t7 * t34 + t5 * t35; 0, 0, 0, 0, 0, 0, 0.2e1 * qJD(3) * t68, (qJ(3) * t68 - t65) * qJD(3), -t23 * t37 - t32 * t34, -t66 + t83, -t25, -t26, 0, t51 * t24 + t40 * t35 - t80, -t51 * t23 - t40 * t34 - t81, t15 * t24 + t3 * t36 + t4 * t30 + t6 * t35 - t80, -t1 * t36 - t19 * t24 + t2 * t37 + t23 * t64 - t8 * t30 + t9 * t32 - t5 * t34 - t7 * t35, t15 * t23 - t3 * t37 - t4 * t32 + t6 * t34 + t81, t1 * t19 + t3 * t15 - t2 * t64 + t6 * t4 + t5 * t9 + t7 * t8; 0, 0, 0, 0, 0, 0, -t82 * qJD(2) ^ 2, t65 * qJD(2), 0, 0, 0, 0, 0, t61, -t12, t61, -t29 - t89, t12, t7 * t30 + (-qJD(5) - t5) * t32 + t67; 0, 0, 0, 0, 0, 0, 0, 0, t87, -t29 + t89, t13, 0, 0, -t40 * t32 + t63, t69 * qJD(4) + t40 * t30 + t62, -t14 * t30 + t63 - t88, pkin(4) * t23 - t24 * qJ(5) + (-t11 + t7) * t32 + (t5 - t76) * t30, t14 * t32 - t6 * t30 + (0.2e1 * qJD(5) - t69) * qJD(4) - t62, -t2 * pkin(4) + t1 * qJ(5) - t5 * t11 - t6 * t14 + t76 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, t13, -qJD(4) ^ 2 - t89, -t7 * qJD(4) + t2 + t88;];
tauc_reg = t16;
