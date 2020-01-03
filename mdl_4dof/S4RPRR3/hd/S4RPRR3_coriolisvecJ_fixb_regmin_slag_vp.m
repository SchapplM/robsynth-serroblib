% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% tauc_reg [4x18]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPRR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:49:18
% EndTime: 2019-12-31 16:49:19
% DurationCPUTime: 0.33s
% Computational Cost: add. (294->68), mult. (749->113), div. (0->0), fcn. (480->6), ass. (0->59)
t34 = sin(pkin(7)) * pkin(1) + pkin(5);
t72 = pkin(6) + t34;
t37 = qJD(3) + qJD(4);
t74 = qJD(4) - t37;
t43 = sin(qJ(3));
t45 = cos(qJ(3));
t53 = t72 * qJD(1);
t14 = t45 * qJD(2) - t53 * t43;
t42 = sin(qJ(4));
t44 = cos(qJ(4));
t25 = t42 * t45 + t44 * t43;
t75 = qJD(1) * t25;
t15 = t43 * qJD(2) + t53 * t45;
t24 = t42 * t43 - t44 * t45;
t8 = t37 * t24;
t73 = t8 * t37;
t64 = qJD(1) * t45;
t57 = t44 * t64;
t65 = qJD(1) * t43;
t58 = t42 * t65;
t18 = -t57 + t58;
t20 = -t42 * t64 - t44 * t65;
t71 = t20 * t18;
t70 = t44 * t15;
t46 = qJD(3) ^ 2;
t69 = t46 * t43;
t68 = t46 * t45;
t67 = t43 ^ 2 - t45 ^ 2;
t66 = qJD(3) * pkin(3);
t35 = -cos(pkin(7)) * pkin(1) - pkin(2);
t28 = qJD(1) * t35;
t63 = t28 * qJD(1);
t61 = qJD(1) * qJD(3);
t60 = pkin(3) * t65;
t59 = t43 * t66;
t13 = t14 + t66;
t56 = -pkin(3) * t37 - t13;
t55 = t45 * t61;
t54 = qJD(3) * t72;
t26 = -t45 * pkin(3) + t35;
t50 = 0.2e1 * qJD(3) * t28;
t11 = t14 * qJD(3);
t12 = t15 * qJD(3);
t21 = t26 * qJD(1);
t49 = -t42 * t11 - t44 * t12 + t21 * t20;
t4 = qJD(4) * t57 - t37 * t58 + t44 * t55;
t48 = t21 * t18 + (t74 * t15 + t12) * t42;
t9 = t37 * t25;
t47 = qJD(1) ^ 2;
t23 = t72 * t45;
t22 = t72 * t43;
t17 = t45 * t54;
t16 = t43 * t54;
t6 = t9 * t37;
t5 = t9 * qJD(1);
t3 = -t18 ^ 2 + t20 ^ 2;
t2 = (-t20 - t75) * t37;
t1 = t18 * t37 + t4;
t7 = [0, 0, 0, 0, 0.2e1 * t43 * t55, -0.2e1 * t67 * t61, t68, -t69, 0, -t34 * t68 + t43 * t50, t34 * t69 + t45 * t50, t20 * t8 + t4 * t25, t8 * t18 + t20 * t9 - t4 * t24 - t25 * t5, -t73, -t6, 0, t26 * t5 + t21 * t9 + (t42 * t16 - t44 * t17 + (t22 * t42 - t23 * t44) * qJD(4)) * t37 + (qJD(1) * t24 + t18) * t59, t26 * t4 - t21 * t8 - (-t44 * t16 - t42 * t17 + (-t22 * t44 - t23 * t42) * qJD(4)) * t37 + (-t20 + t75) * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, -t68, 0, 0, 0, 0, 0, -t6, t73; 0, 0, 0, 0, -t43 * t47 * t45, t67 * t47, 0, 0, 0, -t43 * t63, -t45 * t63, -t71, t3, t1, t2, 0, -t18 * t60 - (-t42 * t14 - t70) * t37 + (t56 * t42 - t70) * qJD(4) + t49, t20 * t60 + (t56 * qJD(4) + t14 * t37 - t11) * t44 + t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, t3, t1, t2, 0, t49 + t74 * (-t42 * t13 - t70), (-t74 * t13 - t11) * t44 + t48;];
tauc_reg = t7;
