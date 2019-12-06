% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tauc_reg [5x18]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRPR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:17:35
% EndTime: 2019-12-05 16:17:37
% DurationCPUTime: 0.41s
% Computational Cost: add. (388->84), mult. (763->151), div. (0->0), fcn. (401->6), ass. (0->73)
t30 = qJD(2) + qJD(3);
t38 = cos(qJ(3));
t66 = pkin(2) * qJD(3);
t55 = qJD(2) * t66;
t17 = t30 * qJD(4) + t38 * t55;
t33 = sin(pkin(9));
t28 = t33 ^ 2;
t34 = cos(pkin(9));
t69 = t34 ^ 2 + t28;
t86 = t69 * t17;
t67 = pkin(2) * qJD(2);
t44 = -t38 * t67 + qJD(4);
t36 = sin(qJ(3));
t85 = t36 * pkin(2);
t84 = t38 * pkin(2);
t19 = t30 * qJ(4) + t36 * t67;
t8 = -t34 * qJD(1) + t33 * t19;
t83 = t8 * t33;
t35 = sin(qJ(5));
t37 = cos(qJ(5));
t20 = -t34 * pkin(4) - t33 * pkin(7) - pkin(3);
t5 = t20 * t30 + t44;
t50 = t36 * t55;
t73 = t34 * t37;
t77 = t28 * t37;
t9 = t33 * qJD(1) + t34 * t19;
t82 = (t35 * t50 + t17 * t73 + (-t35 * t9 + t37 * t5) * qJD(5)) * t34 + t17 * t77;
t24 = t38 * t66 + qJD(4);
t81 = t24 * t30;
t27 = t30 ^ 2;
t80 = t28 * t27;
t79 = t28 * t30;
t78 = t28 * t35;
t76 = t30 * t33;
t75 = t34 * t30;
t74 = t34 * t35;
t72 = t34 * t38;
t23 = -qJD(5) + t75;
t71 = t37 * t23;
t68 = t35 ^ 2 - t37 ^ 2;
t65 = qJD(5) * t33;
t64 = qJD(5) + t23;
t63 = t8 * t76;
t62 = t30 * t77;
t61 = t38 * t79;
t60 = t36 * t66;
t58 = qJD(5) * t79;
t57 = t35 * t65;
t56 = t37 * t65;
t53 = t76 * t85;
t52 = t23 * t57;
t47 = -t35 * t5 - t37 * t9;
t2 = t47 * qJD(5) - t17 * t74 + t37 * t50;
t51 = t17 * t78 - t2 * t34 + t8 * t56;
t49 = t64 * t76;
t48 = t9 * t34 + t83;
t46 = (-qJD(3) + t30) * t67;
t45 = (-qJD(2) - t30) * t66;
t43 = t23 * t34 + t79;
t42 = t36 * t45;
t41 = t36 * t46;
t40 = qJD(4) * t43;
t39 = -t23 ^ 2 - t80;
t25 = qJ(4) + t85;
t21 = t33 * t50;
t18 = -t30 * pkin(3) + t44;
t16 = t57 * t75;
t15 = -0.2e1 * t37 * t35 * t58;
t14 = t20 - t84;
t7 = 0.2e1 * t68 * t58;
t4 = (t23 + t75) * t56;
t3 = t16 + t52;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t23 - t75) * t56, t16 - t52; 0, 0, 0, 0, 0, t42, t38 * t45, t34 * t42, qJD(3) * t53 + t21, t69 * t81 + t86, t48 * t24 + t25 * t86 + (t18 + (-pkin(3) - t84) * qJD(2)) * t60, t15, t7, t3, t4, 0, -(-t24 * t74 + t37 * t60) * t23 + t78 * t81 + (-(-t14 * t35 - t25 * t73) * t23 + t25 * t62) * qJD(5) + t51, (t24 * t73 + t35 * t60) * t23 + t24 * t62 + (t14 * t71 + (-t43 * t25 - t83) * t35) * qJD(5) + t82; 0, 0, 0, 0, 0, t41, t38 * t46, t34 * t41, -qJD(2) * t53 + t21, t44 * t30 * t69 + t86, t48 * qJD(4) + qJ(4) * t86 + (-t48 * t38 + (-pkin(3) * qJD(3) - t18) * t36) * t67, t15, t7, t3, t4, 0, t35 * t40 + (-(-qJ(4) * t73 - t20 * t35) * t23 + qJ(4) * t62) * qJD(5) + ((-t35 * t72 + t36 * t37) * t23 - t35 * t61) * t67 + t51, t37 * t40 + (t20 * t71 + (-t43 * qJ(4) - t83) * t35) * qJD(5) + (-(t35 * t36 + t37 * t72) * t23 - t37 * t61) * t67 + t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69 * t27, -t48 * t30 + t50, 0, 0, 0, 0, 0, t39 * t35, t39 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27 * t35 * t77, -t68 * t80, -t35 * t49, -t37 * t49, 0, t47 * t23 - t37 * t63 + t2, (-t34 * t17 - t64 * t5) * t37 + (t64 * t9 - t50 + t63) * t35;];
tauc_reg = t1;
