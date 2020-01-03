% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RPRR8
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
% Datum: 2019-12-31 16:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPRR8_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR8_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR8_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR8_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:55:15
% EndTime: 2019-12-31 16:55:17
% DurationCPUTime: 0.33s
% Computational Cost: add. (279->77), mult. (659->127), div. (0->0), fcn. (406->4), ass. (0->65)
t58 = 2 * qJD(1);
t37 = qJD(3) + qJD(4);
t77 = qJD(4) - t37;
t44 = (-pkin(1) - pkin(5));
t40 = sin(qJ(4));
t41 = sin(qJ(3));
t42 = cos(qJ(4));
t43 = cos(qJ(3));
t72 = t42 * t43;
t50 = t37 * t72;
t61 = qJD(4) * t40;
t63 = qJD(3) * t41;
t8 = -t40 * t63 - t41 * t61 + t50;
t76 = t8 * t37;
t75 = pkin(6) - t44;
t64 = qJD(1) * t43;
t55 = t42 * t64;
t65 = qJD(1) * t41;
t56 = t40 * t65;
t16 = -t55 + t56;
t21 = t40 * t43 + t42 * t41;
t17 = t21 * qJD(1);
t74 = t16 * t17;
t30 = t44 * qJD(1) + qJD(2);
t14 = -pkin(6) * t65 + t41 * t30;
t73 = t42 * t14;
t45 = qJD(3) ^ 2;
t71 = t45 * t41;
t70 = t45 * t43;
t59 = qJD(1) * qJD(3);
t52 = t41 * t59;
t69 = -qJD(4) * t56 - t40 * t52;
t68 = t41 ^ 2 - t43 ^ 2;
t46 = qJD(1) ^ 2;
t67 = -t45 - t46;
t66 = t46 * qJ(2);
t62 = qJD(3) * t43;
t60 = qJ(2) * qJD(3);
t57 = pkin(3) * t64;
t54 = qJD(2) * t58;
t26 = t75 * t43;
t15 = -pkin(6) * t64 + t43 * t30;
t11 = qJD(3) * pkin(3) + t15;
t53 = -pkin(3) * t37 - t11;
t35 = t41 * pkin(3) + qJ(2);
t51 = pkin(6) * qJD(1) - t30;
t31 = pkin(3) * t62 + qJD(2);
t12 = t51 * t63;
t13 = t51 * t62;
t27 = t35 * qJD(1);
t48 = t42 * t12 + t40 * t13 + t27 * t16;
t47 = t14 * t61 + (-t14 * t37 - t12) * t40 + t27 * t17;
t7 = t37 * t21;
t4 = t7 * qJD(1);
t25 = t75 * t41;
t24 = t31 * qJD(1);
t22 = -t40 * t41 + t72;
t20 = qJD(3) * t26;
t19 = t75 * t63;
t6 = t7 * t37;
t5 = qJD(1) * t50 + t69;
t3 = t16 ^ 2 - t17 ^ 2;
t2 = -t69 + (-t16 - t55) * t37;
t1 = t17 * t37 - t4;
t9 = [0, 0, 0, 0, t54, qJ(2) * t54, -0.2e1 * t43 * t52, 0.2e1 * t68 * t59, -t71, -t70, 0, -t44 * t71 + (qJD(2) * t41 + t43 * t60) * t58, -t44 * t70 + (qJD(2) * t43 - t41 * t60) * t58, t16 * t7 - t4 * t22, t16 * t8 + t7 * t17 + t4 * t21 - t22 * t5, -t6, -t76, 0, t31 * t17 + t35 * t5 + t24 * t21 + t27 * t8 + (t42 * t19 + t40 * t20 + (t25 * t42 + t26 * t40) * qJD(4)) * t37, -t31 * t16 - t35 * t4 + t24 * t22 - t27 * t7 - (t40 * t19 - t42 * t20 + (t25 * t40 - t26 * t42) * qJD(4)) * t37; 0, 0, 0, 0, -t46, -t66, 0, 0, 0, 0, 0, t67 * t41, t67 * t43, 0, 0, 0, 0, 0, -qJD(1) * t17 - t6, qJD(1) * t16 - t76; 0, 0, 0, 0, 0, 0, t43 * t46 * t41, -t68 * t46, 0, 0, 0, -t43 * t66, t41 * t66, -t74, t3, t1, t2, 0, -t17 * t57 - (-t40 * t15 - t73) * t37 + (t53 * t40 - t73) * qJD(4) + t48, t16 * t57 + (t53 * qJD(4) + t15 * t37 + t13) * t42 + t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, t3, t1, t2, 0, t48 + t77 * (-t40 * t11 - t73), (-t77 * t11 + t13) * t42 + t47;];
tauc_reg = t9;
