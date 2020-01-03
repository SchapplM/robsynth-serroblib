% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% tauc_reg [4x21]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRPR8_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR8_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:08:23
% EndTime: 2019-12-31 17:08:25
% DurationCPUTime: 0.63s
% Computational Cost: add. (352->107), mult. (915->176), div. (0->0), fcn. (520->4), ass. (0->83)
t48 = sin(qJ(4));
t49 = sin(qJ(2));
t50 = cos(qJ(4));
t51 = cos(qJ(2));
t21 = -t51 * t48 + t49 * t50;
t78 = qJD(2) * t51;
t102 = qJD(4) * t21 + t48 * t78;
t20 = t48 * t49 + t50 * t51;
t59 = t20 * qJD(4);
t75 = qJD(1) * qJD(2);
t70 = t51 * t75;
t71 = t49 * t75;
t1 = -qJD(1) * t59 + t48 * t71 + t50 * t70;
t13 = t20 * qJD(1);
t43 = qJD(2) - qJD(4);
t101 = -t13 * t43 + t1;
t100 = -0.2e1 * t75;
t80 = qJD(1) * t51;
t81 = qJD(1) * t49;
t15 = -t48 * t80 + t50 * t81;
t2 = qJD(1) * t102 - t50 * t71;
t99 = t15 * t43 + t2;
t97 = qJD(4) + t43;
t96 = pkin(2) + pkin(3);
t95 = pkin(5) - pkin(6);
t93 = t15 * t13;
t54 = qJD(1) ^ 2;
t89 = t51 * t54;
t53 = qJD(2) ^ 2;
t88 = t53 * t49;
t87 = t53 * t51;
t46 = t49 ^ 2;
t84 = -t51 ^ 2 + t46;
t83 = qJ(3) * t51;
t82 = qJD(2) * pkin(2);
t62 = pkin(2) * t49 - t83;
t77 = t49 * qJD(3);
t12 = qJD(2) * t62 - t77;
t7 = qJD(1) * t12;
t69 = t49 * qJ(3) + pkin(1);
t29 = -pkin(2) * t51 - t69;
t17 = qJD(1) * t29;
t79 = qJD(2) * t49;
t41 = pkin(5) * t81;
t76 = -pkin(6) * t81 + qJD(3) + t41;
t45 = qJD(2) * qJ(3);
t74 = t49 * t89;
t42 = pkin(5) * t80;
t32 = t95 * t51;
t72 = t13 ^ 2 - t15 ^ 2;
t68 = 0.2e1 * t17;
t67 = pkin(1) * t100;
t66 = qJD(3) - t82;
t65 = t43 ^ 2;
t24 = t95 * t79;
t25 = -pkin(6) * t80 + t42;
t16 = t25 + t45;
t9 = -qJD(2) * t96 + t76;
t64 = -t50 * t16 - t48 * t9;
t63 = t48 * t16 - t50 * t9;
t61 = -pkin(5) * t53 - 0.2e1 * t7;
t60 = -t49 * t96 + t83;
t19 = t51 * t96 + t69;
t44 = qJD(2) * qJD(3);
t10 = -qJD(1) * t24 + t44;
t40 = pkin(5) * t70;
t18 = -pkin(6) * t70 + t40;
t8 = t19 * qJD(1);
t58 = -t50 * t10 + t8 * t13 - t48 * t18;
t57 = -t48 * t10 - t8 * t15 + t50 * t18;
t6 = qJD(2) * t60 + t77;
t27 = -pkin(5) * t71 + t44;
t28 = t41 + t66;
t30 = t42 + t45;
t55 = t27 * t51 + (t28 * t51 + (-t30 + t42) * t49) * qJD(2);
t31 = t95 * t49;
t26 = qJD(2) * t32;
t22 = t62 * qJD(1);
t11 = t60 * qJD(1);
t5 = t6 * qJD(1);
t4 = qJD(2) * t20 - t59;
t3 = -t50 * t79 + t102;
t14 = [0, 0, 0, 0.2e1 * t49 * t70, t84 * t100, t87, -t88, 0, -pkin(5) * t87 + t49 * t67, pkin(5) * t88 + t51 * t67, t51 * t61 + t68 * t79, t55, t49 * t61 - t68 * t78, pkin(5) * t55 + t17 * t12 + t7 * t29, t1 * t21 + t15 * t4, -t1 * t20 - t13 * t4 - t15 * t3 - t2 * t21, -t4 * t43, t3 * t43, 0, t6 * t13 + t19 * t2 + t5 * t20 + t8 * t3 - (t48 * t24 + t50 * t26 + (-t31 * t48 - t32 * t50) * qJD(4)) * t43, t6 * t15 + t19 * t1 + t5 * t21 + t8 * t4 + (-t50 * t24 + t48 * t26 + (t31 * t50 - t32 * t48) * qJD(4)) * t43; 0, 0, 0, -t74, t84 * t54, 0, 0, 0, t54 * pkin(1) * t49, pkin(1) * t89, (-t17 * t49 + t22 * t51) * qJD(1), ((t30 - t45) * t49 + (-t28 + t66) * t51) * qJD(1), 0.2e1 * t44 + (t17 * t51 + t22 * t49) * qJD(1), t27 * qJ(3) + t30 * qJD(3) - t17 * t22 + (t30 * t49 + (-t28 - t82) * t51) * qJD(1) * pkin(5), -t93, t72, -t101, t99, 0, -t11 * t13 + (t50 * t25 + t76 * t48) * t43 + (-(-qJ(3) * t50 + t48 * t96) * t43 - t64) * qJD(4) - t57, -t11 * t15 + (-t48 * t25 + t76 * t50) * t43 + ((-qJ(3) * t48 - t50 * t96) * t43 - t63) * qJD(4) - t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, 0, -t46 * t54 - t53, -qJD(2) * t30 + t17 * t81 + t40, 0, 0, 0, 0, 0, -t13 * t81 - t48 * t65, -t15 * t81 - t50 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, -t72, t101, -t99, 0, t64 * t97 + t57, t63 * t97 + t58;];
tauc_reg = t14;
