% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tauc_reg [5x17]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:34
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:34:20
% EndTime: 2022-01-20 10:34:22
% DurationCPUTime: 0.33s
% Computational Cost: add. (682->84), mult. (1469->121), div. (0->0), fcn. (862->8), ass. (0->70)
t54 = qJD(5) ^ 2;
t47 = cos(pkin(9));
t38 = pkin(2) * t47 + pkin(3);
t49 = sin(qJ(4));
t52 = cos(qJ(4));
t46 = sin(pkin(9));
t92 = pkin(2) * t46;
t61 = t38 * t49 + t52 * t92;
t53 = cos(qJ(2));
t50 = sin(qJ(2));
t86 = t47 * t50;
t58 = pkin(1) * (-t46 * t53 - t86);
t29 = qJD(1) * t58;
t87 = t46 * t50;
t57 = pkin(1) * (t47 * t53 - t87);
t31 = qJD(1) * t57;
t43 = qJD(1) + qJD(2);
t41 = qJD(4) + t43;
t74 = (-t61 * qJD(4) - t29 * t52 + t31 * t49) * t41;
t94 = (pkin(8) + t61) * t54 - t74;
t79 = pkin(1) * qJD(1);
t36 = pkin(2) * t43 + t53 * t79;
t76 = t50 * t79;
t18 = t47 * t36 - t46 * t76;
t16 = pkin(3) * t43 + t18;
t30 = qJD(2) * t58;
t22 = qJD(1) * t30;
t32 = qJD(2) * t57;
t23 = qJD(1) * t32;
t19 = t36 * t46 + t47 * t76;
t84 = t49 * t19;
t56 = -(qJD(4) * t16 + t23) * t52 + qJD(4) * t84 - t49 * t22;
t11 = t16 * t49 + t19 * t52;
t3 = qJD(4) * t11 - t52 * t22 + t49 * t23;
t48 = sin(qJ(5));
t51 = cos(qJ(5));
t10 = t16 * t52 - t84;
t91 = t41 * pkin(4);
t8 = -t10 - t91;
t78 = qJD(5) * t8;
t93 = t3 * t48 + t51 * t78;
t39 = pkin(1) * t53 + pkin(2);
t70 = -pkin(1) * t87 + t47 * t39;
t26 = pkin(3) + t70;
t33 = pkin(1) * t86 + t39 * t46;
t64 = t26 * t49 + t33 * t52;
t90 = (qJD(4) * t64 - t52 * t30 + t49 * t32) * t41;
t89 = t11 * t41;
t85 = t48 * t51;
t83 = t54 * t48;
t60 = t38 * t52 - t49 * t92;
t81 = -t60 * qJD(4) + t29 * t49 + t31 * t52;
t80 = t48 ^ 2 - t51 ^ 2;
t77 = 0.2e1 * qJD(5) * t41;
t75 = -t41 * t8 + t56;
t69 = (-qJD(2) + t43) * t79;
t68 = pkin(1) * qJD(2) * (-qJD(1) - t43);
t67 = pkin(8) * t54 - t89;
t66 = (pkin(8) + t64) * t54 + t90;
t65 = t26 * t52 - t33 * t49;
t63 = qJD(5) * (t10 - t91);
t4 = qJD(4) * t65 + t49 * t30 + t52 * t32;
t62 = qJD(5) * ((-pkin(4) - t65) * t41 - t4);
t59 = qJD(5) * ((-pkin(4) - t60) * t41 + t81);
t42 = t54 * t51;
t40 = t41 ^ 2;
t35 = t77 * t85;
t21 = t80 * t77;
t6 = t48 * t78;
t1 = [0, 0, 0, 0, t50 * t68, t53 * t68, t18 * t30 + t19 * t32 + t22 * t70 + t23 * t33, 0, -t3 - t90, -t4 * t41 + t56, t35, -t21, t42, -t83, 0, t6 + t48 * t62 + (-t3 - t66) * t51, t48 * t66 + t51 * t62 + t93; 0, 0, 0, 0, t50 * t69, t53 * t69, -t18 * t29 - t19 * t31 + (t22 * t47 + t23 * t46) * pkin(2), 0, -t3 + t74, t41 * t81 + t56, t35, -t21, t42, -t83, 0, t6 + t48 * t59 + (-t3 - t94) * t51, t48 * t94 + t51 * t59 + t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83, -t42; 0, 0, 0, 0, 0, 0, 0, 0, -t3 + t89, t10 * t41 + t56, t35, -t21, t42, -t83, 0, t6 + t48 * t63 + (-t3 - t67) * t51, t48 * t67 + t51 * t63 + t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40 * t85, t80 * t40, 0, 0, 0, t75 * t48, t75 * t51;];
tauc_reg = t1;
