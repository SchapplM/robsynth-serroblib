% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% tauc_reg [4x16]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRRR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:22:14
% EndTime: 2019-12-31 17:22:15
% DurationCPUTime: 0.24s
% Computational Cost: add. (314->69), mult. (665->98), div. (0->0), fcn. (340->6), ass. (0->60)
t30 = sin(qJ(3));
t34 = cos(qJ(2));
t31 = sin(qJ(2));
t33 = cos(qJ(3));
t67 = t31 * t33;
t45 = t30 * t34 + t67;
t63 = pkin(1) * qJD(1);
t14 = t45 * t63;
t26 = qJD(1) + qJD(2);
t24 = qJD(3) + t26;
t35 = qJD(4) ^ 2;
t61 = qJD(3) * t30;
t76 = (pkin(2) * t61 - t14) * t24 + (t30 * pkin(2) + pkin(7)) * t35;
t57 = t34 * t63;
t17 = t26 * pkin(2) + t57;
t52 = qJD(2) * t57;
t58 = t31 * t63;
t53 = t30 * t58;
t65 = (qJD(2) + qJD(3)) * t53;
t39 = -(qJD(3) * t17 + t52) * t33 + t65;
t60 = qJD(3) * t33;
t75 = t45 * qJD(2) + t31 * t60;
t29 = sin(qJ(4));
t38 = t75 * pkin(1);
t56 = t17 * t61;
t3 = qJD(1) * t38 + t56;
t32 = cos(qJ(4));
t73 = t24 * pkin(3);
t9 = t33 * t17 - t53;
t8 = -t9 - t73;
t62 = qJD(4) * t8;
t74 = t3 * t29 + t32 * t62;
t22 = t34 * pkin(1) + pkin(2);
t72 = (t22 * t61 + t38) * t24;
t71 = (t30 * t17 + t33 * t58) * t24;
t69 = t29 * t32;
t68 = t30 * t31;
t66 = t35 * t29;
t64 = t29 ^ 2 - t32 ^ 2;
t59 = 0.2e1 * qJD(4) * t24;
t54 = -t24 * t8 + t39;
t50 = (-qJD(2) + t26) * t63;
t49 = pkin(1) * qJD(2) * (-qJD(1) - t26);
t48 = pkin(7) * t35 - t71;
t47 = (pkin(1) * t67 + t30 * t22 + pkin(7)) * t35 + t72;
t46 = qJD(4) * (t9 - t73);
t44 = t33 * t34 - t68;
t4 = t22 * t60 + (t44 * qJD(2) - t31 * t61) * pkin(1);
t43 = qJD(4) * ((pkin(1) * t68 - t33 * t22 - pkin(3)) * t24 - t4);
t42 = (-pkin(2) * t24 - t17) * qJD(3);
t15 = t44 * t63;
t40 = qJD(4) * (-pkin(2) * t60 + (-t33 * pkin(2) - pkin(3)) * t24 + t15);
t37 = t75 * t63;
t36 = -t37 - t56;
t25 = t35 * t32;
t23 = t24 ^ 2;
t16 = t59 * t69;
t11 = t64 * t59;
t6 = t29 * t62;
t1 = [0, 0, 0, 0, t31 * t49, t34 * t49, 0, t36 - t72, -t4 * t24 + t39, t16, -t11, t25, -t66, 0, t6 + t29 * t43 + (-t3 - t47) * t32, t47 * t29 + t32 * t43 + t74; 0, 0, 0, 0, t31 * t50, t34 * t50, 0, t14 * t24 + t30 * t42 - t37, t15 * t24 + (t42 - t52) * t33 + t65, t16, -t11, t25, -t66, 0, t6 + t29 * t40 + (-t3 - t76) * t32, t76 * t29 + t32 * t40 + t74; 0, 0, 0, 0, 0, 0, 0, t36 + t71, t9 * t24 + t39, t16, -t11, t25, -t66, 0, t6 + t29 * t46 + (-t3 - t48) * t32, t48 * t29 + t32 * t46 + t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23 * t69, t64 * t23, 0, 0, 0, t54 * t29, t54 * t32;];
tauc_reg = t1;
