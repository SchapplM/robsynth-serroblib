% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% tauc_reg [5x16]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRRP1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:40:12
% EndTime: 2019-12-05 16:40:15
% DurationCPUTime: 0.38s
% Computational Cost: add. (424->88), mult. (769->132), div. (0->0), fcn. (370->4), ass. (0->72)
t80 = qJ(5) + pkin(7);
t36 = sin(qJ(4));
t38 = cos(qJ(4));
t33 = qJD(2) + qJD(3);
t37 = sin(qJ(3));
t70 = pkin(2) * qJD(2);
t58 = t37 * t70;
t50 = t80 * t33 + t58;
t7 = t38 * qJD(1) - t50 * t36;
t8 = t36 * qJD(1) + t50 * t38;
t68 = qJD(4) * pkin(4);
t4 = t7 + t68;
t79 = t4 - t7;
t39 = cos(qJ(3));
t78 = t39 * pkin(2);
t77 = t33 * t36;
t76 = t37 * t38;
t40 = qJD(4) ^ 2;
t75 = t40 * t36;
t30 = t40 * t38;
t57 = t39 * t70;
t20 = -t33 * pkin(3) - t57;
t69 = pkin(2) * qJD(3);
t51 = qJD(2) * t69;
t25 = t37 * t51;
t64 = t38 * qJD(4);
t73 = t20 * t64 + t36 * t25;
t65 = t36 * qJD(4);
t56 = t33 * t65;
t13 = pkin(4) * t56 + t25;
t34 = t36 ^ 2;
t35 = t38 ^ 2;
t72 = -t34 - t35;
t71 = t34 - t35;
t26 = t37 * pkin(2) + pkin(7);
t67 = -qJ(5) - t26;
t63 = -qJD(2) - t33;
t62 = -qJD(3) + t33;
t61 = t39 * t69;
t60 = t37 * t69;
t59 = pkin(4) * t65;
t55 = t33 * t64;
t54 = t39 * t65;
t53 = -t38 * pkin(4) - pkin(3);
t47 = t39 * t51;
t41 = qJD(5) * t33 + t47;
t2 = t7 * qJD(4) + t41 * t38;
t3 = -qJD(4) * t8 - t41 * t36;
t52 = t2 * t38 - t3 * t36;
t49 = qJD(4) * t80;
t48 = qJD(4) * t67;
t46 = -t36 * t8 - t38 * t4;
t45 = t36 * t4 - t38 * t8;
t43 = -t37 * t77 + t39 * t64;
t42 = -t20 * t33 - t47;
t32 = t33 ^ 2;
t31 = t38 * qJ(5);
t28 = t38 * qJD(5);
t27 = -pkin(3) - t78;
t23 = t38 * pkin(7) + t31;
t22 = t80 * t36;
t18 = 0.2e1 * t36 * t55;
t17 = t38 * t26 + t31;
t16 = t67 * t36;
t14 = t20 * t65;
t12 = -t36 * qJD(5) - t38 * t49;
t11 = -t36 * t49 + t28;
t10 = -0.2e1 * t71 * t33 * qJD(4);
t9 = t53 * t33 + qJD(5) - t57;
t6 = (-qJD(5) - t61) * t36 + t38 * t48;
t5 = t36 * t48 + t38 * t61 + t28;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75, -t30, 0, -t45 * qJD(4) + t2 * t36 + t3 * t38; 0, 0, 0, 0, 0, -t33 * t60 - t25, t63 * t61, t18, t10, t30, -t75, 0, t27 * t56 - t26 * t30 + t14 + (t63 * t76 - t54) * t69, t26 * t75 + t27 * t55 - t43 * t69 + t73, (-t36 * t6 + t38 * t5) * t33 + ((-t16 * t38 - t17 * t36) * t33 + t46) * qJD(4) + t52, t2 * t17 + t8 * t5 + t3 * t16 + t4 * t6 + t13 * (t53 - t78) + t9 * (t59 + t60); 0, 0, 0, 0, 0, t33 * t58 - t25, t62 * t57, t18, t10, t30, -t75, 0, -pkin(3) * t56 - pkin(7) * t30 + t14 + (t62 * t76 + t54) * t70, -pkin(3) * t55 + pkin(7) * t75 + t43 * t70 + t73, t46 * qJD(4) + (t11 * t38 - t12 * t36 + (t22 * t38 - t23 * t36) * qJD(4) + t72 * t57) * t33 + t52, t2 * t23 + t8 * t11 - t3 * t22 + t4 * t12 + t13 * t53 + t9 * t59 + (-t37 * t9 + t45 * t39) * t70; 0, 0, 0, 0, 0, 0, 0, -t36 * t32 * t38, t71 * t32, 0, 0, 0, t42 * t36, t42 * t38, (-t68 + t79) * t38 * t33, t79 * t8 + (-t9 * t77 + t3) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72 * t32, t45 * t33 + t13;];
tauc_reg = t1;
