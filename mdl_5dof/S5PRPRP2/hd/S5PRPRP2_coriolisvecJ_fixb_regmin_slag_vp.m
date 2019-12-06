% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRPRP2
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
% tauc_reg [5x17]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRPRP2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:31:04
% EndTime: 2019-12-05 15:31:07
% DurationCPUTime: 0.42s
% Computational Cost: add. (390->90), mult. (992->156), div. (0->0), fcn. (599->4), ass. (0->73)
t32 = cos(pkin(8));
t33 = sin(qJ(4));
t34 = cos(qJ(4));
t31 = sin(pkin(8));
t19 = -t32 * pkin(3) - t31 * pkin(6) - pkin(2);
t69 = qJ(5) * t31;
t47 = -t19 + t69;
t81 = -t34 * t32 * qJ(3) + t47 * t33;
t58 = qJ(3) * qJD(2);
t17 = t31 * qJD(1) + t32 * t58;
t64 = qJD(4) * t33;
t14 = t19 * qJD(2) + qJD(3);
t65 = qJD(3) * t34;
t23 = t32 * t65;
t63 = qJD(4) * t34;
t74 = -qJD(2) * t23 - t14 * t63;
t38 = -t17 * t64 - t74;
t57 = qJ(5) * qJD(4);
t62 = qJD(5) * t33;
t67 = qJD(2) * t31;
t1 = (-t34 * t57 - t62) * t67 + t38;
t40 = -t33 * t14 - t34 * t17;
t37 = t40 * qJD(4);
t66 = qJD(3) * t33;
t52 = t32 * t66;
t61 = qJD(5) * t34;
t2 = t37 + (-t52 + (t33 * t57 - t61) * t31) * qJD(2);
t60 = t32 * qJD(2);
t24 = -qJD(4) + t60;
t12 = t34 * t14;
t50 = qJ(5) * t67;
t6 = -t33 * t17 - t34 * t50 + t12;
t3 = -t24 * pkin(4) + t6;
t7 = -t33 * t50 - t40;
t42 = t3 * t33 - t34 * t7;
t80 = -t42 * qJD(4) + t1 * t33 + t2 * t34;
t28 = t32 ^ 2;
t27 = t31 ^ 2;
t79 = 0.2e1 * t27;
t78 = t3 - t6;
t26 = t32 * qJD(1);
t16 = t31 * t58 - t26;
t77 = t16 * t31;
t76 = t24 * t32;
t35 = qJD(2) ^ 2;
t75 = t27 * t35;
t73 = t19 * t63 + t23;
t56 = qJD(2) * qJD(4);
t49 = t34 * t56;
t13 = t31 * pkin(4) * t49 + qJD(3) * t67;
t72 = t27 + t28;
t29 = t33 ^ 2;
t30 = t34 ^ 2;
t71 = t29 - t30;
t70 = qJ(3) * t33;
t68 = qJD(2) * t27;
t59 = qJD(4) + t24;
t54 = t31 * t64;
t53 = t31 * t63;
t51 = qJ(3) * t64;
t48 = pkin(4) * t33 + qJ(3);
t46 = t72 * qJD(2);
t45 = t24 * t54;
t43 = t3 * t34 + t33 * t7;
t41 = t59 * t67;
t39 = t17 * t32 + t77;
t36 = -t24 ^ 2 - t75;
t18 = t54 * t60;
t10 = t48 * t67 + qJD(5) - t26;
t8 = -t47 * t34 + (-pkin(4) - t70) * t32;
t5 = t81 * qJD(4) - t31 * t61 - t52;
t4 = -t31 * t62 + (-t32 * t70 - t34 * t69) * qJD(4) + t73;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t24 - t60) * t53, t18 - t45, 0, -t13 * t32 + (-t43 * qJD(4) + t1 * t34 - t2 * t33) * t31; 0, 0, 0, 0, 0, 0, 0.2e1 * qJD(3) * t46, (qJ(3) * t46 + t39) * qJD(3), -0.2e1 * t27 * t33 * t49, t71 * t56 * t79, t18 + t45, (t24 + t60) * t53, 0, (t76 + (t79 + t28) * qJD(2)) * t66 + ((t14 * t32 + t19 * t24) * t33 + ((t68 + t76) * qJ(3) + t39) * t34) * qJD(4), (-t32 * t51 + t73) * t24 + t38 * t32 - t16 * t54 + (-t51 + 0.2e1 * t65) * t68, ((-t33 * t4 - t34 * t5 + (t33 * t8 + t34 * t81) * qJD(4)) * qJD(2) - t80) * t31, -t1 * t81 + t2 * t8 + t3 * t5 + t7 * t4 + (t13 * t48 + t10 * (pkin(4) * t63 + qJD(3))) * t31; 0, 0, 0, 0, 0, 0, -t72 * t35, -t39 * qJD(2), 0, 0, 0, 0, 0, t36 * t33, t36 * t34, 0, (-t10 * t31 + t42 * t32) * qJD(2) + t80; 0, 0, 0, 0, 0, 0, 0, 0, t34 * t33 * t75, -t71 * t75, -t33 * t41, -t34 * t41, 0, t40 * t24 + t37 + (-t34 * t77 - t52) * qJD(2), -t12 * t24 + (t16 * t67 + t59 * t17) * t33 + t74, (pkin(4) * qJD(4) - t78) * t33 * t67, t78 * t7 + (-t10 * t34 * t67 + t2) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t29 - t30) * t75, t43 * t67 + t13;];
tauc_reg = t9;
