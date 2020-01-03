% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% 
% Output:
% tauc_reg [4x20]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPRR9_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR9_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:56:28
% EndTime: 2019-12-31 16:56:30
% DurationCPUTime: 0.51s
% Computational Cost: add. (348->108), mult. (816->187), div. (0->0), fcn. (468->4), ass. (0->70)
t24 = sin(qJ(3));
t52 = t24 * qJD(1);
t18 = qJD(4) + t52;
t75 = qJD(4) - t18;
t47 = 0.2e1 * qJD(1);
t25 = cos(qJ(4));
t26 = cos(qJ(3));
t56 = qJD(1) * t26;
t43 = t25 * t56;
t23 = sin(qJ(4));
t53 = t23 * qJD(3);
t13 = t43 + t53;
t49 = qJD(1) * qJD(3);
t68 = t23 * t24;
t4 = qJD(4) * t13 - t49 * t68;
t27 = -pkin(1) - pkin(5);
t54 = qJD(4) * t23;
t45 = t26 * t54;
t51 = t25 * qJD(3);
t30 = -t24 * t51 - t45;
t3 = t30 * qJD(1) + qJD(4) * t51;
t74 = t3 * t23;
t17 = t27 * qJD(1) + qJD(2);
t9 = -qJD(3) * pkin(3) - t26 * t17;
t73 = t9 * t23;
t72 = t9 * t25;
t11 = t23 * t56 - t51;
t71 = t11 * t18;
t70 = t13 * t18;
t69 = t23 * t18;
t67 = t23 * t26;
t66 = t24 * t27;
t65 = t25 * t18;
t64 = t25 * t26;
t63 = t26 * t27;
t62 = t27 * t18;
t28 = qJD(3) ^ 2;
t61 = t28 * t24;
t60 = t28 * t26;
t22 = t26 ^ 2;
t59 = t24 ^ 2 - t22;
t29 = qJD(1) ^ 2;
t58 = -t28 - t29;
t57 = t29 * qJ(2);
t55 = qJD(3) * t26;
t50 = qJ(2) * qJD(3);
t48 = t18 * t64;
t46 = t17 * t55;
t44 = qJD(4) * t65;
t42 = qJD(2) * t47;
t41 = t26 * t49;
t40 = t18 + t52;
t39 = t11 + t51;
t38 = -t13 + t53;
t37 = qJD(4) * t24 + qJD(1);
t36 = pkin(3) * t26 + pkin(6) * t24;
t15 = t24 * pkin(3) - t26 * pkin(6) + qJ(2);
t7 = t15 * qJD(1);
t8 = qJD(3) * pkin(6) + t24 * t17;
t35 = t23 * t8 - t25 * t7;
t2 = t23 * t7 + t25 * t8;
t34 = qJD(1) * t22 - t18 * t24;
t33 = -pkin(6) * t55 + t9 * t24;
t32 = t23 * t15 + t25 * t66;
t31 = t25 * t15 - t23 * t66;
t10 = t36 * qJD(3) + qJD(2);
t14 = t36 * qJD(1);
t6 = t10 * qJD(1);
t5 = t25 * t6;
t1 = [0, 0, 0, 0, t42, qJ(2) * t42, -0.2e1 * t24 * t41, 0.2e1 * t59 * t49, -t61, -t60, 0, -t27 * t61 + (qJD(2) * t24 + t26 * t50) * t47, -t27 * t60 + (qJD(2) * t26 - t24 * t50) * t47, t30 * t13 + t3 * t64, (t11 * t25 + t13 * t23) * t24 * qJD(3) + (-t74 - t25 * t4 + (t11 * t23 - t13 * t25) * qJD(4)) * t26, -t18 * t45 + t3 * t24 + (t13 * t26 + t34 * t25) * qJD(3), -t26 * t44 - t4 * t24 + (-t11 * t26 - t34 * t23) * qJD(3), t40 * t55, t10 * t65 - t4 * t63 + t5 * t24 + (-t32 * t18 - t2 * t24 + t9 * t64) * qJD(4) + ((t27 * t11 - t73) * t24 + (t31 * qJD(1) - t23 * t62 - t35) * t26) * qJD(3), -t3 * t63 + (-t10 * t18 - t24 * t6) * t23 + (-t31 * t18 + t35 * t24 - t9 * t67) * qJD(4) + ((t27 * t13 - t72) * t24 + (-t32 * qJD(1) - t25 * t62 - t2) * t26) * qJD(3); 0, 0, 0, 0, -t29, -t57, 0, 0, 0, 0, 0, t58 * t24, t58 * t26, 0, 0, 0, 0, 0, -t26 * t4 - t37 * t65 + (t24 * t11 - t40 * t67) * qJD(3), -t26 * t3 + t37 * t69 + (-t48 + (t13 - t43) * t24) * qJD(3); 0, 0, 0, 0, 0, 0, t26 * t29 * t24, -t59 * t29, 0, 0, 0, -t26 * t57, t24 * t57, t13 * t65 + t74, (t3 - t71) * t25 + (-t4 - t70) * t23, t44 + (t24 * t65 + t38 * t26) * qJD(1), -t18 * t54 + (-t18 * t68 + t39 * t26) * qJD(1), -t18 * t56, -t14 * t65 - pkin(3) * t4 + (t18 * t67 - t39 * t24) * t17 + (-pkin(6) * t65 + t73) * qJD(4) + (t33 * t23 + t26 * t35) * qJD(1), t14 * t69 - pkin(3) * t3 + (t38 * t24 + t48) * t17 + (pkin(6) * t69 + t72) * qJD(4) + (t2 * t26 + t33 * t25) * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13 * t11, -t11 ^ 2 + t13 ^ 2, t3 + t71, -t4 + t70, t41, -t9 * t13 - t75 * t2 - t23 * t46 + t5, t9 * t11 - t23 * t6 - t25 * t46 + t75 * t35;];
tauc_reg = t1;
