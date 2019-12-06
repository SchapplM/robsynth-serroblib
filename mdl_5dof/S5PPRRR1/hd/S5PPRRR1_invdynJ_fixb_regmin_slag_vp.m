% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% tau_reg [5x15]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PPRRR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:12:52
% EndTime: 2019-12-05 15:12:54
% DurationCPUTime: 0.47s
% Computational Cost: add. (580->100), mult. (1137->138), div. (0->0), fcn. (960->14), ass. (0->77)
t50 = pkin(9) + qJ(3);
t47 = qJ(4) + t50;
t41 = sin(t47);
t55 = sin(pkin(8));
t57 = cos(pkin(8));
t79 = g(1) * t57 + g(2) * t55;
t115 = t79 * t41;
t54 = sin(pkin(9));
t56 = cos(pkin(9));
t60 = sin(qJ(3));
t63 = cos(qJ(3));
t28 = t54 * t63 + t56 * t60;
t26 = t28 * qJD(3);
t109 = -t54 * t60 + t56 * t63;
t78 = t109 * qJDD(1);
t15 = qJDD(3) * pkin(3) - qJD(1) * t26 + t78;
t106 = t28 * qJDD(1);
t25 = t109 * qJD(3);
t16 = qJD(1) * t25 + t106;
t23 = t109 * qJD(1);
t21 = qJD(3) * pkin(3) + t23;
t59 = sin(qJ(4));
t62 = cos(qJ(4));
t24 = t28 * qJD(1);
t42 = cos(t47);
t89 = qJD(4) * t59;
t82 = g(3) * t41 + t24 * t89 + t79 * t42;
t114 = -(qJD(4) * t21 + t16) * t62 - t59 * t15 + t82;
t113 = -g(3) * t42 + t15 * t62 - t59 * t16;
t49 = qJDD(3) + qJDD(4);
t100 = t49 * pkin(4);
t92 = t62 * t24;
t11 = t21 * t59 + t92;
t107 = t11 * qJD(4);
t111 = -t100 + t107 - t113;
t13 = t23 * t59 + t92;
t43 = pkin(3) * t59 + pkin(7);
t98 = t62 * pkin(3);
t44 = -pkin(4) - t98;
t51 = qJD(3) + qJD(4);
t64 = qJD(5) ^ 2;
t105 = t43 * t64 + t44 * t49 + (pkin(3) * t89 - t13) * t51;
t99 = t51 * pkin(4);
t97 = t11 * t51;
t95 = t59 * t24;
t61 = cos(qJ(5));
t93 = t61 * t49;
t58 = sin(qJ(5));
t52 = t58 ^ 2;
t90 = -t61 ^ 2 + t52;
t88 = t61 * qJD(5);
t10 = t21 * t62 - t95;
t8 = -t10 - t99;
t87 = t111 * t58 + t8 * t88;
t86 = t8 * qJD(5) * t58 + t115 * t61;
t84 = -pkin(3) * t51 - t21;
t18 = t109 * t59 + t28 * t62;
t76 = t109 * t62 - t28 * t59;
t77 = t76 * t49 - (qJD(4) * t18 + t59 * t25 + t62 * t26) * t51;
t73 = g(1) * t55 - g(2) * t57 - qJDD(2);
t71 = pkin(7) * t64 - t100 - t97;
t70 = t18 * t64 - t77;
t69 = t115 + t113;
t68 = -pkin(7) * qJDD(5) + (t10 - t99) * qJD(5);
t4 = qJD(4) * t76 + t62 * t25 - t59 * t26;
t67 = -qJDD(5) * t18 + (-t51 * t76 - t4) * qJD(5);
t66 = -t49 * pkin(7) - t8 * t51 + t114;
t14 = t23 * t62 - t95;
t65 = -qJDD(5) * t43 + (-qJD(4) * t98 + t44 * t51 + t14) * qJD(5);
t48 = t51 ^ 2;
t46 = cos(t50);
t45 = sin(t50);
t36 = qJDD(5) * t61 - t58 * t64;
t35 = qJDD(5) * t58 + t61 * t64;
t22 = 0.2e1 * t51 * t58 * t88 + t49 * t52;
t19 = -0.2e1 * qJD(5) * t51 * t90 + 0.2e1 * t58 * t93;
t1 = [qJDD(1) - g(3), -g(3) + (t54 ^ 2 + t56 ^ 2) * qJDD(1), 0, -qJD(3) * t26 + qJDD(3) * t109, -qJD(3) * t25 - qJDD(3) * t28, 0, t77, -t18 * t49 - t4 * t51, 0, 0, 0, 0, 0, t58 * t67 - t61 * t70, t58 * t70 + t61 * t67; 0, -t73, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t35; 0, 0, qJDD(3), -g(3) * t46 + t45 * t79 + t78, g(3) * t45 + t79 * t46 - t106, t49, t49 * t98 + t13 * t51 + (t59 * t84 - t92) * qJD(4) + t69, t14 * t51 + (-pkin(3) * t49 - t15) * t59 + (qJD(4) * t84 - t16) * t62 + t82, t22, t19, t35, t36, 0, t65 * t58 + (-t105 - t111) * t61 + t86, t65 * t61 + (t105 - t115) * t58 + t87; 0, 0, 0, 0, 0, t49, t69 + t97 - t107, t10 * t51 + t114, t22, t19, t35, t36, 0, t68 * t58 + (-t71 - t111) * t61 + t86, t68 * t61 + (-t115 + t71) * t58 + t87; 0, 0, 0, 0, 0, 0, 0, 0, -t58 * t48 * t61, t90 * t48, t58 * t49, t93, qJDD(5), t58 * t66 - t61 * t73, t58 * t73 + t61 * t66;];
tau_reg = t1;
