% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RRPP5
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
% Datum: 2019-12-31 17:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRPP5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:00:37
% EndTime: 2019-12-31 17:00:38
% DurationCPUTime: 0.33s
% Computational Cost: add. (279->97), mult. (708->138), div. (0->0), fcn. (294->2), ass. (0->76)
t42 = cos(qJ(2));
t64 = qJD(1) * t42;
t32 = pkin(5) * t64;
t18 = pkin(3) * t64 + t32;
t78 = -qJD(4) - t18;
t60 = qJD(1) * qJD(2);
t77 = -0.2e1 * t60;
t61 = qJD(2) * qJ(3);
t12 = t61 - t78;
t40 = -pkin(2) - qJ(4);
t76 = qJD(2) * t40;
t75 = pkin(3) + pkin(5);
t41 = sin(qJ(2));
t65 = qJD(1) * t41;
t31 = pkin(5) * t65;
t67 = qJD(2) * pkin(2);
t21 = qJD(3) + t31 - t67;
t74 = t21 * t42;
t23 = -t32 - t61;
t73 = t23 * t41;
t44 = qJD(1) ^ 2;
t72 = t42 * t44;
t43 = qJD(2) ^ 2;
t71 = t43 * t41;
t70 = t43 * t42;
t56 = t42 * t60;
t27 = pkin(5) * t56;
t69 = pkin(3) * t56 + t27;
t38 = t41 ^ 2;
t39 = t42 ^ 2;
t68 = t38 - t39;
t66 = t41 * qJ(3);
t14 = t40 * t42 - pkin(1) - t66;
t4 = qJD(1) * t14;
t51 = -t42 * pkin(2) - t66;
t22 = -pkin(1) + t51;
t13 = qJD(1) * t22;
t63 = t41 * qJD(3);
t16 = -pkin(3) * t65 - t31;
t62 = qJD(3) - t16;
t59 = t41 * t72;
t25 = t75 * t42;
t58 = t75 * qJD(2);
t57 = t41 * t60;
t28 = pkin(2) * t57;
t50 = -qJ(3) * t42 + qJ(4) * t41;
t46 = t50 * qJD(2) - t42 * qJD(4) - t63;
t1 = t46 * qJD(1) + t28;
t34 = t41 * t67;
t2 = t34 + t46;
t55 = -qJD(1) * t2 - t1;
t54 = 0.2e1 * t4;
t53 = pkin(1) * t77;
t17 = t41 * t58;
t49 = -0.2e1 * qJD(2) * t13;
t47 = -t42 * t61 - t63;
t11 = t34 + t47;
t5 = t47 * qJD(1) + t28;
t48 = pkin(5) * t43 + qJD(1) * t11 + t5;
t37 = qJD(2) * qJD(3);
t20 = pkin(5) * t57 - t37;
t45 = -t20 * t42 + (t74 + (t23 + t32) * t41) * qJD(2);
t36 = 0.2e1 * t37;
t35 = pkin(2) * t65;
t30 = qJD(3) * t64;
t26 = -t38 * t44 - t43;
t24 = t75 * t41;
t19 = qJD(2) * t25;
t15 = -qJ(3) * t64 + t35;
t10 = t50 * qJD(1) + t35;
t9 = -qJD(2) * qJD(4) + t69;
t8 = -qJD(1) * t17 + t37;
t7 = t62 + t76;
t6 = t13 * t65;
t3 = t4 * t64;
t29 = [0, 0, 0, 0.2e1 * t41 * t56, t68 * t77, t70, -t71, 0, -pkin(5) * t70 + t41 * t53, pkin(5) * t71 + t42 * t53, t45, t41 * t49 + t48 * t42, -t48 * t41 + t42 * t49, t45 * pkin(5) + t13 * t11 + t5 * t22, t9 * t41 + t8 * t42 + (-t12 * t41 + t42 * t7) * qJD(2) + (-t17 * t42 + t19 * t41 + (t24 * t42 - t25 * t41) * qJD(2)) * qJD(1), t55 * t41 + (-t54 * t42 - t17) * qJD(2), t55 * t42 + (t54 * t41 - t19) * qJD(2), t1 * t14 - t12 * t17 + t7 * t19 + t4 * t2 + t9 * t24 + t8 * t25; 0, 0, 0, -t59, t68 * t44, 0, 0, 0, t44 * pkin(1) * t41, pkin(1) * t72, t30 + (t51 * qJD(2) - t73 - t74) * qJD(1), -t15 * t64 + t6, t36 + (t13 * t42 + t15 * t41) * qJD(1), -t20 * qJ(3) - t23 * qJD(3) - t13 * t15 + (-t73 + (-t21 - t67) * t42) * qJD(1) * pkin(5), t30 + (-t16 - t7 + t76) * t64, -t16 * qJD(2) + t3 + t36 + (t10 - t58) * t65, (0.2e1 * qJD(4) + t18) * qJD(2) + (t10 * t42 - t4 * t41) * qJD(1) - t69, t8 * qJ(3) - t4 * t10 + t62 * t12 + t9 * t40 + t78 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, t26, t23 * qJD(2) + t27 + t6, 0, t26, -t59, t4 * t65 + (-qJD(4) - t12) * qJD(2) + t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, -t39 * t44 - t43, t3 + t37 + (-t75 * t65 + t7) * qJD(2);];
tauc_reg = t29;
