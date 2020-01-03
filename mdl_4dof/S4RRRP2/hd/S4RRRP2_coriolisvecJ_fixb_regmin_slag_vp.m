% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% tauc_reg [4x15]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRRP2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:13:09
% EndTime: 2019-12-31 17:13:11
% DurationCPUTime: 0.33s
% Computational Cost: add. (347->82), mult. (648->129), div. (0->0), fcn. (301->4), ass. (0->74)
t79 = qJ(4) + pkin(6);
t67 = qJD(3) * pkin(3);
t35 = sin(qJ(3));
t32 = qJD(1) + qJD(2);
t36 = sin(qJ(2));
t69 = pkin(1) * qJD(1);
t57 = t36 * t69;
t49 = t79 * t32 + t57;
t7 = t49 * t35;
t6 = -t7 + t67;
t78 = t6 + t7;
t38 = cos(qJ(2));
t77 = t38 * pkin(1);
t76 = t32 * t35;
t37 = cos(qJ(3));
t75 = t36 * t37;
t39 = qJD(3) ^ 2;
t74 = t39 * t35;
t29 = t39 * t37;
t56 = t38 * t69;
t20 = -t32 * pkin(2) - t56;
t68 = pkin(1) * qJD(2);
t50 = qJD(1) * t68;
t25 = t36 * t50;
t64 = qJD(3) * t37;
t72 = t20 * t64 + t35 * t25;
t65 = qJD(3) * t35;
t55 = t32 * t65;
t13 = pkin(3) * t55 + t25;
t33 = t35 ^ 2;
t34 = t37 ^ 2;
t71 = -t33 - t34;
t70 = t33 - t34;
t26 = t36 * pkin(1) + pkin(6);
t66 = -qJ(4) - t26;
t63 = qJD(3) * t38;
t62 = -qJD(1) - t32;
t61 = -qJD(2) + t32;
t60 = t38 * t68;
t59 = t36 * t68;
t58 = pkin(3) * t65;
t54 = t32 * t64;
t53 = t35 * t63;
t52 = -t37 * pkin(3) - pkin(2);
t46 = t38 * t50;
t40 = qJD(4) * t32 + t46;
t43 = qJD(3) * t49;
t2 = -t35 * t43 + t40 * t37;
t3 = -t40 * t35 - t37 * t43;
t51 = t2 * t37 - t3 * t35;
t48 = qJD(3) * t79;
t47 = qJD(3) * t66;
t8 = t49 * t37;
t45 = -t35 * t8 - t37 * t6;
t44 = t35 * t6 - t37 * t8;
t42 = -t36 * t76 + t37 * t63;
t41 = -t20 * t32 - t46;
t31 = t32 ^ 2;
t30 = t37 * qJ(4);
t28 = t37 * qJD(4);
t27 = -pkin(2) - t77;
t23 = t37 * pkin(6) + t30;
t22 = t79 * t35;
t18 = 0.2e1 * t35 * t54;
t17 = t37 * t26 + t30;
t16 = t66 * t35;
t14 = t20 * t65;
t12 = -t35 * qJD(4) - t37 * t48;
t11 = -t35 * t48 + t28;
t10 = -0.2e1 * t70 * t32 * qJD(3);
t9 = t52 * t32 + qJD(4) - t56;
t5 = (-qJD(4) - t60) * t35 + t37 * t47;
t4 = t35 * t47 + t37 * t60 + t28;
t1 = [0, 0, 0, 0, -t32 * t59 - t25, t62 * t60, t18, t10, t29, -t74, 0, t27 * t55 - t26 * t29 + t14 + (t62 * t75 - t53) * t68, t26 * t74 + t27 * t54 - t42 * t68 + t72, (-t35 * t5 + t37 * t4) * t32 + ((-t16 * t37 - t17 * t35) * t32 + t45) * qJD(3) + t51, t2 * t17 + t8 * t4 + t3 * t16 + t6 * t5 + t13 * (t52 - t77) + t9 * (t58 + t59); 0, 0, 0, 0, t32 * t57 - t25, t61 * t56, t18, t10, t29, -t74, 0, -pkin(2) * t55 - pkin(6) * t29 + t14 + (t61 * t75 + t53) * t69, -pkin(2) * t54 + pkin(6) * t74 + t42 * t69 + t72, t45 * qJD(3) + (t11 * t37 - t12 * t35 + (t22 * t37 - t23 * t35) * qJD(3) + t71 * t56) * t32 + t51, t2 * t23 + t8 * t11 - t3 * t22 + t6 * t12 + t13 * t52 + t9 * t58 + (-t36 * t9 + t44 * t38) * t69; 0, 0, 0, 0, 0, 0, -t35 * t31 * t37, t70 * t31, 0, 0, 0, t41 * t35, t41 * t37, (-t67 + t78) * t37 * t32, t78 * t8 + (-t9 * t76 + t3) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71 * t31, t44 * t32 + t13;];
tauc_reg = t1;
