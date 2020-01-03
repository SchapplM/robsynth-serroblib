% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% tau_reg [4x14]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRPR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:01:35
% EndTime: 2019-12-31 17:01:36
% DurationCPUTime: 0.34s
% Computational Cost: add. (346->88), mult. (621->132), div. (0->0), fcn. (360->12), ass. (0->73)
t50 = qJ(1) + qJ(2);
t42 = pkin(7) + t50;
t34 = cos(t42);
t57 = cos(qJ(2));
t88 = t57 * pkin(1);
t41 = qJDD(1) * t88;
t46 = qJDD(1) + qJDD(2);
t54 = sin(qJ(2));
t92 = pkin(1) * t54;
t74 = qJD(1) * t92;
t14 = t46 * pkin(2) - qJD(2) * t74 + t41;
t51 = sin(pkin(7));
t52 = cos(pkin(7));
t79 = qJD(1) * t57;
t71 = qJD(2) * t79;
t76 = qJDD(1) * t54;
t4 = -(t71 + t76) * pkin(1) * t51 + t52 * t14;
t95 = -t46 * pkin(3) + g(2) * t34 - t4;
t43 = sin(t50);
t44 = cos(t50);
t94 = g(1) * t43 - g(2) * t44;
t33 = sin(t42);
t91 = g(1) * t33;
t53 = sin(qJ(4));
t56 = cos(qJ(4));
t47 = qJD(1) + qJD(2);
t73 = pkin(1) * t79;
t23 = t47 * pkin(2) + t73;
t11 = t52 * t23 - t51 * t74;
t9 = -t47 * pkin(3) - t11;
t87 = t9 * qJD(4) * t53 + t56 * t91;
t86 = t51 * t54;
t84 = t52 * t54;
t83 = t56 * t46;
t32 = pkin(1) * t84;
t40 = pkin(2) + t88;
t82 = t51 * t40 + t32;
t81 = g(1) * t44 + g(2) * t43;
t48 = t53 ^ 2;
t80 = -t56 ^ 2 + t48;
t78 = t56 * qJD(4);
t77 = qJDD(3) - g(3);
t75 = t95 * t53 + t9 * t78;
t5 = t52 * pkin(1) * t71 + qJDD(1) * t32 + t51 * t14;
t70 = qJD(1) * (-qJD(2) + t47);
t69 = qJD(2) * (-qJD(1) - t47);
t68 = t41 + t94;
t67 = -pkin(1) * t86 + t52 * t40;
t66 = pkin(1) * (t52 * t57 - t86);
t16 = -pkin(3) - t67;
t17 = pkin(6) + t82;
t19 = (t51 * t57 + t84) * qJD(2) * pkin(1);
t59 = qJD(4) ^ 2;
t64 = t16 * t46 + t17 * t59 + t19 * t47;
t31 = t52 * t74;
t18 = t51 * t73 + t31;
t35 = t51 * pkin(2) + pkin(6);
t36 = -t52 * pkin(2) - pkin(3);
t63 = -t18 * t47 + t35 * t59 + t36 * t46;
t62 = -t46 * pkin(6) + g(1) * t34 + g(2) * t33 - t9 * t47 - t5;
t21 = qJD(2) * t66;
t61 = -qJDD(4) * t17 + (t16 * t47 - t21) * qJD(4);
t20 = qJD(1) * t66;
t60 = -qJDD(4) * t35 + (t36 * t47 + t20) * qJD(4);
t58 = cos(qJ(1));
t55 = sin(qJ(1));
t45 = t47 ^ 2;
t25 = qJDD(4) * t56 - t59 * t53;
t24 = qJDD(4) * t53 + t59 * t56;
t15 = 0.2e1 * t53 * t47 * t78 + t48 * t46;
t12 = t51 * t23 + t31;
t8 = -0.2e1 * t80 * t47 * qJD(4) + 0.2e1 * t53 * t83;
t1 = [qJDD(1), g(1) * t55 - g(2) * t58, g(1) * t58 + g(2) * t55, t46, (t46 * t57 + t54 * t69) * pkin(1) + t68, ((-qJDD(1) - t46) * t54 + t57 * t69) * pkin(1) + t81, t5 * t82 + t12 * t21 + t4 * t67 - t11 * t19 - g(1) * (-t55 * pkin(1) - pkin(2) * t43) - g(2) * (t58 * pkin(1) + pkin(2) * t44), t15, t8, t24, t25, 0, t61 * t53 + (-t64 - t95) * t56 + t87, t61 * t56 + (t64 - t91) * t53 + t75; 0, 0, 0, t46, t70 * t92 + t68, (t57 * t70 - t76) * pkin(1) + t81, t11 * t18 - t12 * t20 + (t4 * t52 + t5 * t51 + t94) * pkin(2), t15, t8, t24, t25, 0, t60 * t53 + (-t63 - t95) * t56 + t87, t60 * t56 + (t63 - t91) * t53 + t75; 0, 0, 0, 0, 0, 0, t77, 0, 0, 0, 0, 0, t25, -t24; 0, 0, 0, 0, 0, 0, 0, -t53 * t45 * t56, t80 * t45, t53 * t46, t83, qJDD(4), t62 * t53 + t77 * t56, -t77 * t53 + t62 * t56;];
tau_reg = t1;
