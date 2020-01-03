% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RRRP2
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% tau_reg [4x15]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRRP2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP2_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:13:11
% EndTime: 2019-12-31 17:13:12
% DurationCPUTime: 0.52s
% Computational Cost: add. (574->133), mult. (869->177), div. (0->0), fcn. (439->8), ass. (0->90)
t107 = qJ(4) + pkin(6);
t45 = qJDD(1) + qJDD(2);
t100 = t45 * pkin(2);
t49 = qJ(1) + qJ(2);
t42 = cos(t49);
t102 = g(2) * t42;
t52 = sin(qJ(2));
t90 = pkin(1) * qJD(1);
t78 = t52 * t90;
t55 = cos(qJ(2));
t98 = t55 * pkin(1);
t93 = -qJD(2) * t78 + qJDD(1) * t98;
t106 = -t93 - t100 + t102;
t41 = sin(t49);
t105 = g(1) * t42 + g(2) * t41;
t34 = g(1) * t41;
t104 = t102 - t34;
t89 = qJD(3) * pkin(3);
t51 = sin(qJ(3));
t46 = qJD(1) + qJD(2);
t72 = t107 * t46 + t78;
t9 = t72 * t51;
t8 = -t9 + t89;
t103 = t8 + t9;
t54 = cos(qJ(3));
t101 = g(3) * t54;
t99 = t46 * pkin(2);
t97 = t51 * t45;
t96 = t54 * t45;
t87 = qJD(1) * t55;
t77 = pkin(1) * t87;
t23 = -t77 - t99;
t84 = qJD(3) * t51;
t95 = t23 * t84 + t54 * t34;
t47 = t51 ^ 2;
t48 = t54 ^ 2;
t92 = -t47 - t48;
t91 = t47 - t48;
t36 = t52 * pkin(1) + pkin(6);
t88 = -qJ(4) - t36;
t86 = qJD(2) * t52;
t85 = qJD(2) * t55;
t83 = qJDD(1) * t52;
t14 = t45 * pkin(6) + (qJD(1) * t85 + t83) * pkin(1);
t63 = qJ(4) * t45 + qJD(4) * t46 + t14;
t64 = qJD(3) * t72;
t3 = -t51 * t64 + t63 * t54;
t82 = t3 * t54 - t105;
t81 = t23 * qJD(3) * t54 + t106 * t51;
t80 = pkin(1) * t85;
t79 = pkin(3) * t84;
t76 = t46 * t86;
t75 = t46 * t84;
t37 = t54 * pkin(3) + pkin(2);
t73 = t107 * t41 + t42 * t37;
t71 = qJD(3) * t107;
t70 = qJD(3) * t88;
t69 = t46 * t78;
t10 = t72 * t54;
t67 = -t10 * t54 + t51 * t8;
t66 = t107 * t42 - t41 * t37;
t65 = -t93 + t104;
t62 = -t23 * t46 + t105 - t14;
t57 = qJD(3) ^ 2;
t61 = pkin(6) * t57 - t100 - t69;
t4 = pkin(3) * t75 - t37 * t45 + qJDD(4) - t93;
t38 = -pkin(2) - t98;
t60 = pkin(1) * t76 + t36 * t57 + t38 * t45;
t59 = -pkin(6) * qJDD(3) + (t77 - t99) * qJD(3);
t58 = -qJDD(3) * t36 + (t38 * t46 - t80) * qJD(3);
t56 = cos(qJ(1));
t53 = sin(qJ(1));
t44 = t46 ^ 2;
t43 = t54 * qJ(4);
t40 = t54 * qJD(4);
t28 = t54 * pkin(6) + t43;
t27 = t107 * t51;
t26 = qJDD(3) * t54 - t57 * t51;
t25 = qJDD(3) * t51 + t57 * t54;
t21 = t54 * t36 + t43;
t20 = t88 * t51;
t17 = -t51 * qJD(4) - t54 * t71;
t16 = -t51 * t71 + t40;
t15 = t47 * t45 + 0.2e1 * t54 * t75;
t12 = -t37 * t46 + qJD(4) - t77;
t7 = -0.2e1 * t91 * t46 * qJD(3) + 0.2e1 * t51 * t96;
t6 = (-qJD(4) - t80) * t51 + t54 * t70;
t5 = t51 * t70 + t54 * t80 + t40;
t2 = qJDD(3) * pkin(3) - t63 * t51 - t54 * t64;
t1 = [qJDD(1), g(1) * t53 - g(2) * t56, g(1) * t56 + g(2) * t53, t45, (t45 * t55 - t76) * pkin(1) - t65, ((-qJDD(1) - t45) * t52 + (-qJD(1) - t46) * t85) * pkin(1) + t105, t15, t7, t25, t26, 0, t58 * t51 + (-t60 - t106) * t54 + t95, t58 * t54 + (t60 - t34) * t51 + t81, (t21 * t45 + t46 * t5 + (-t20 * t46 - t8) * qJD(3)) * t54 + (-t20 * t45 - t46 * t6 - t2 + (-t21 * t46 - t10) * qJD(3)) * t51 + t82, t3 * t21 + t10 * t5 + t2 * t20 + t8 * t6 + t4 * (-t37 - t98) + t12 * (pkin(1) * t86 + t79) - g(1) * (-t53 * pkin(1) + t66) - g(2) * (t56 * pkin(1) + t73); 0, 0, 0, t45, -t65 + t69, (-t83 + (-qJD(2) + t46) * t87) * pkin(1) + t105, t15, t7, t25, t26, 0, t59 * t51 + (-t61 - t106) * t54 + t95, t59 * t54 + (t61 - t34) * t51 + t81, (-qJD(3) * t8 + t28 * t45) * t54 + (-qJD(3) * t10 + t27 * t45 - t2) * t51 + (t16 * t54 - t17 * t51 + (t27 * t54 - t28 * t51) * qJD(3) + t92 * t77) * t46 + t82, t3 * t28 + t10 * t16 - t2 * t27 + t8 * t17 - t4 * t37 + t12 * t79 - g(1) * t66 - g(2) * t73 + (-t12 * t52 + t67 * t55) * t90; 0, 0, 0, 0, 0, 0, -t51 * t44 * t54, t91 * t44, t97, t96, qJDD(3), t62 * t51 - t101, g(3) * t51 + t62 * t54, -pkin(3) * t97 + (-t89 + t103) * t54 * t46, t103 * t10 + (-t101 + t2 + (-t12 * t46 + t105) * t51) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92 * t44, t67 * t46 + t104 + t4;];
tau_reg = t1;
