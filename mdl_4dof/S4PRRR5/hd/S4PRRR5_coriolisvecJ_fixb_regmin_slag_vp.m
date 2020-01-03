% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4PRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% tauc_reg [4x14]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PRRR5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:33:46
% EndTime: 2019-12-31 16:33:47
% DurationCPUTime: 0.24s
% Computational Cost: add. (188->52), mult. (431->83), div. (0->0), fcn. (284->6), ass. (0->52)
t69 = 2 * qJD(4);
t32 = cos(qJ(2));
t51 = t32 * qJD(1);
t19 = qJD(2) * pkin(2) + t51;
t29 = sin(qJ(2));
t28 = sin(qJ(3));
t31 = cos(qJ(3));
t15 = t28 * t32 + t31 * t29;
t38 = t15 * qJD(2);
t52 = qJD(3) * t31;
t35 = (t29 * t52 + t38) * qJD(1);
t53 = qJD(3) * t28;
t3 = t19 * t53 + t35;
t24 = qJD(2) + qJD(3);
t14 = t28 * t29 - t31 * t32;
t68 = t14 * t24;
t12 = t15 * qJD(1);
t33 = qJD(4) ^ 2;
t66 = (pkin(2) * t53 - t12) * t24 + (t28 * pkin(2) + pkin(6)) * t33;
t45 = qJD(2) * t51;
t54 = qJD(1) * t29;
t47 = t28 * t54;
t57 = t24 * t47;
t65 = -(qJD(3) * t19 + t45) * t31 + t57;
t27 = sin(qJ(4));
t30 = cos(qJ(4));
t63 = t24 * pkin(3);
t9 = t31 * t19 - t47;
t8 = -t9 - t63;
t55 = qJD(4) * t8;
t64 = t3 * t27 + t30 * t55;
t62 = (t15 * qJD(3) + t38) * t24;
t61 = (t28 * t19 + t31 * t54) * t24;
t59 = t27 * t30;
t58 = t33 * t27;
t56 = t27 ^ 2 - t30 ^ 2;
t50 = t24 * t69;
t46 = -t24 * t8 + t65;
t43 = pkin(6) * t33 - t61;
t42 = t15 * t33 + t62;
t41 = qJD(4) * (t9 - t63);
t40 = t68 * t69;
t39 = (-pkin(2) * t24 - t19) * qJD(3);
t13 = t14 * qJD(1);
t36 = qJD(4) * (-pkin(2) * t52 + (-t31 * pkin(2) - pkin(3)) * t24 - t13);
t34 = qJD(2) ^ 2;
t23 = t24 ^ 2;
t22 = t33 * t30;
t16 = t50 * t59;
t11 = t56 * t50;
t6 = t27 * t55;
t1 = [0, 0, -t34 * t29, -t34 * t32, 0, -t62, t68 * t24, 0, 0, 0, 0, 0, t27 * t40 - t42 * t30, t42 * t27 + t30 * t40; 0, 0, 0, 0, 0, t12 * t24 + t28 * t39 - t35, -t13 * t24 + (t39 - t45) * t31 + t57, t16, -t11, t22, -t58, 0, t6 + t27 * t36 + (-t3 - t66) * t30, t66 * t27 + t30 * t36 + t64; 0, 0, 0, 0, 0, t61 - t3, t9 * t24 + t65, t16, -t11, t22, -t58, 0, t6 + t27 * t41 + (-t3 - t43) * t30, t43 * t27 + t30 * t41 + t64; 0, 0, 0, 0, 0, 0, 0, -t23 * t59, t56 * t23, 0, 0, 0, t46 * t27, t46 * t30;];
tauc_reg = t1;
