% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% tauc_reg [4x18]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRPP4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:59:18
% EndTime: 2019-12-31 16:59:19
% DurationCPUTime: 0.35s
% Computational Cost: add. (278->98), mult. (716->146), div. (0->0), fcn. (303->2), ass. (0->76)
t77 = pkin(2) + pkin(3);
t54 = t77 * qJD(2);
t36 = sin(qJ(2));
t65 = qJD(1) * t36;
t29 = pkin(5) * t65;
t59 = qJ(4) * qJD(1);
t15 = t36 * t59 - t29;
t60 = qJD(3) - t15;
t7 = -t54 + t60;
t37 = cos(qJ(2));
t51 = t36 * qJ(3) + pkin(1);
t14 = t77 * t37 + t51;
t66 = qJD(1) * t14;
t4 = qJD(4) + t66;
t67 = qJD(4) + t4;
t78 = t67 * t36;
t70 = pkin(5) * qJD(1);
t30 = t37 * t70;
t17 = -t37 * t59 + t30;
t33 = qJD(2) * qJ(3);
t12 = t17 + t33;
t40 = qJD(1) ^ 2;
t76 = t37 * t40;
t39 = qJD(2) ^ 2;
t75 = t39 * t36;
t74 = t39 * t37;
t73 = pkin(5) - qJ(4);
t57 = qJD(1) * qJD(2);
t52 = t37 * t57;
t63 = t36 * qJD(3);
t72 = qJ(3) * t52 + qJD(1) * t63;
t34 = t36 ^ 2;
t35 = t37 ^ 2;
t71 = t34 - t35;
t69 = qJ(3) * t37;
t68 = qJD(2) * pkin(2);
t20 = -t37 * pkin(2) - t51;
t13 = qJD(1) * t20;
t64 = qJD(2) * t36;
t62 = t36 * qJD(4);
t61 = t37 * qJD(4);
t58 = qJ(4) * qJD(2);
t56 = t36 * t76;
t55 = t37 * t58;
t23 = t73 * t37;
t53 = t36 * t57;
t1 = -t77 * t53 + t72;
t42 = -t77 * t36 + t69;
t3 = t42 * qJD(2) + t63;
t50 = qJD(1) * t3 + t1;
t49 = t4 + t66;
t48 = 0.2e1 * t13;
t47 = qJD(3) - t68;
t46 = -0.2e1 * t53;
t45 = 0.2e1 * t52;
t44 = pkin(2) * t36 - t69;
t11 = t44 * qJD(2) - t63;
t6 = pkin(2) * t53 - t72;
t43 = -pkin(5) * t39 - qJD(1) * t11 - t6;
t32 = qJD(2) * qJD(3);
t18 = -pkin(5) * t53 + t32;
t19 = t29 + t47;
t21 = t30 + t33;
t41 = t18 * t37 + (t19 * t37 + (-t21 + t30) * t36) * qJD(2);
t31 = 0.2e1 * t32;
t27 = pkin(5) * t52;
t25 = qJ(4) * t53;
t24 = -t34 * t40 - t39;
t22 = t73 * t36;
t16 = t44 * qJD(1);
t10 = qJD(2) * t23 - t62;
t9 = -t73 * t64 - t61;
t8 = t42 * qJD(1);
t5 = t27 + (-t55 - t62) * qJD(1);
t2 = t25 + t32 + (-pkin(5) * t64 - t61) * qJD(1);
t26 = [0, 0, 0, t36 * t45, -0.2e1 * t71 * t57, t74, -t75, 0, pkin(1) * t46 - pkin(5) * t74, -0.2e1 * pkin(1) * t52 + pkin(5) * t75, t43 * t37 + t48 * t64, t41, -t48 * t37 * qJD(2) + t43 * t36, t41 * pkin(5) + t13 * t11 + t6 * t20, t50 * t37 + (-t49 * t36 - t10) * qJD(2), t50 * t36 + (t49 * t37 + t9) * qJD(2), -t2 * t37 - t5 * t36 + (t12 * t36 - t37 * t7) * qJD(2) + (-t10 * t36 - t37 * t9 + (-t22 * t37 + t23 * t36) * qJD(2)) * qJD(1), t1 * t14 + t7 * t10 + t12 * t9 + t2 * t23 + t5 * t22 + t4 * t3; 0, 0, 0, -t56, t71 * t40, 0, 0, 0, t40 * pkin(1) * t36, pkin(1) * t76, (-t13 * t36 + t16 * t37) * qJD(1), ((t21 - t33) * t36 + (-t19 + t47) * t37) * qJD(1), t31 + (t13 * t37 + t16 * t36) * qJD(1), t18 * qJ(3) + t21 * qJD(3) - t13 * t16 + (t21 * t36 + (-t19 - t68) * t37) * t70, t17 * qJD(2) - t27 + ((-t8 + t58) * t37 + t78) * qJD(1), -t15 * qJD(2) + t25 + t31 + (-t67 * t37 + (-pkin(5) * qJD(2) - t8) * t36) * qJD(1), 0, t2 * qJ(3) + t60 * t12 - t7 * t17 - t4 * t8 - t5 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, 0, t24, -t21 * qJD(2) + t13 * t65 + t27, -t56, t24, 0, -t12 * qJD(2) + t27 + (-t55 - t78) * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, t45, (-t34 - t35) * t40, (t12 * t37 + (t7 - t54) * t36) * qJD(1) + t72;];
tauc_reg = t26;
