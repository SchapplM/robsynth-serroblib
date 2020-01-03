% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% tauc_reg [5x17]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPRP1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:25:49
% EndTime: 2020-01-03 11:25:51
% DurationCPUTime: 0.42s
% Computational Cost: add. (496->93), mult. (1179->161), div. (0->0), fcn. (705->6), ass. (0->74)
t26 = sin(pkin(7)) * pkin(1) + qJ(3);
t35 = cos(pkin(8));
t37 = sin(qJ(4));
t38 = cos(qJ(4));
t33 = sin(pkin(8));
t18 = -cos(pkin(7)) * pkin(1) - pkin(2) - t33 * pkin(6) - t35 * pkin(3);
t72 = qJ(5) * t33;
t51 = -t18 + t72;
t84 = -t38 * t35 * t26 + t51 * t37;
t20 = t26 * qJD(1);
t16 = t33 * qJD(2) + t35 * t20;
t66 = qJD(4) * t37;
t13 = t18 * qJD(1) + qJD(3);
t67 = qJD(3) * t38;
t24 = t35 * t67;
t65 = qJD(4) * t38;
t76 = -qJD(1) * t24 - t13 * t65;
t42 = -t16 * t66 - t76;
t60 = qJ(5) * qJD(4);
t64 = qJD(5) * t37;
t70 = qJD(1) * t33;
t1 = (-t38 * t60 - t64) * t70 + t42;
t44 = -t37 * t13 - t38 * t16;
t41 = t44 * qJD(4);
t68 = qJD(3) * t37;
t54 = t35 * t68;
t63 = qJD(5) * t38;
t2 = t41 + (-t54 + (t37 * t60 - t63) * t33) * qJD(1);
t62 = t35 * qJD(1);
t25 = -qJD(4) + t62;
t12 = t38 * t13;
t53 = qJ(5) * t70;
t6 = -t37 * t16 - t38 * t53 + t12;
t3 = -t25 * pkin(4) + t6;
t7 = -t37 * t53 - t44;
t46 = t3 * t37 - t38 * t7;
t83 = -t46 * qJD(4) + t1 * t37 + t2 * t38;
t30 = t35 ^ 2;
t29 = t33 ^ 2;
t82 = 0.2e1 * t29;
t81 = t3 - t6;
t28 = t35 * qJD(2);
t15 = t33 * t20 - t28;
t80 = t15 * t33;
t79 = t25 * t35;
t78 = t26 * t37;
t39 = qJD(1) ^ 2;
t77 = t29 * t39;
t75 = t18 * t65 + t24;
t59 = qJD(1) * qJD(4);
t52 = t38 * t59;
t17 = t33 * pkin(4) * t52 + qJD(3) * t70;
t74 = t29 + t30;
t31 = t37 ^ 2;
t32 = t38 ^ 2;
t73 = t31 - t32;
t71 = qJD(1) * t29;
t69 = qJD(1) * t37;
t61 = qJD(4) + t25;
t57 = t33 * t66;
t56 = t33 * t65;
t55 = t35 * t66;
t50 = qJD(1) * t74;
t49 = t25 * t57;
t47 = t3 * t38 + t37 * t7;
t45 = t61 * t70;
t43 = t16 * t35 + t80;
t40 = -t25 ^ 2 - t77;
t19 = t55 * t70;
t10 = qJD(5) - t28 + (pkin(4) * t69 + t20) * t33;
t8 = -t51 * t38 + (-pkin(4) - t78) * t35;
t5 = t84 * qJD(4) - t33 * t63 - t54;
t4 = -t33 * t64 + (-t35 * t78 - t38 * t72) * qJD(4) + t75;
t9 = [0, 0, 0, 0, 0, 0, 0.2e1 * qJD(3) * t50, (t26 * t50 + t43) * qJD(3), -0.2e1 * t29 * t37 * t52, t73 * t59 * t82, t19 + t49, (t25 + t62) * t56, 0, (t79 + (t82 + t30) * qJD(1)) * t68 + ((t13 * t35 + t18 * t25) * t37 + ((t71 + t79) * t26 + t43) * t38) * qJD(4), (-t26 * t55 + t75) * t25 + t42 * t35 - t15 * t57 + (-t26 * t66 + 0.2e1 * t67) * t71, ((-t37 * t4 - t38 * t5 + (t37 * t8 + t38 * t84) * qJD(4)) * qJD(1) - t83) * t33, -t1 * t84 + t2 * t8 + t3 * t5 + t7 * t4 + (t17 * (pkin(4) * t37 + t26) + t10 * (pkin(4) * t65 + qJD(3))) * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t25 - t62) * t56, t19 - t49, 0, -t17 * t35 + (-t47 * qJD(4) + t1 * t38 - t2 * t37) * t33; 0, 0, 0, 0, 0, 0, -t74 * t39, -t43 * qJD(1), 0, 0, 0, 0, 0, t40 * t37, t40 * t38, 0, (-t10 * t33 + t46 * t35) * qJD(1) + t83; 0, 0, 0, 0, 0, 0, 0, 0, t38 * t37 * t77, -t73 * t77, -t37 * t45, -t38 * t45, 0, t44 * t25 + t41 + (-t38 * t80 - t54) * qJD(1), -t12 * t25 + (t15 * t70 + t61 * t16) * t37 + t76, (pkin(4) * qJD(4) - t81) * t33 * t69, t81 * t7 + (-t10 * t38 * t70 + t2) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t31 - t32) * t77, t47 * t70 + t17;];
tauc_reg = t9;
