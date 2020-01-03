% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4PRRR6
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
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4PRRR6_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR6_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR6_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR6_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:35:06
% EndTime: 2019-12-31 16:35:08
% DurationCPUTime: 0.60s
% Computational Cost: add. (1275->148), mult. (2636->203), div. (0->0), fcn. (1769->8), ass. (0->96)
t86 = sin(qJ(3));
t103 = qJD(2) * t86;
t85 = sin(qJ(4));
t88 = cos(qJ(4));
t89 = cos(qJ(3));
t51 = -t88 * t89 * qJD(2) + t85 * t103;
t53 = (t89 * t85 + t86 * t88) * qJD(2);
t37 = t53 * t51;
t77 = qJDD(3) + qJDD(4);
t115 = -t37 + t77;
t118 = t115 * t85;
t117 = t115 * t88;
t101 = t86 * qJDD(2);
t99 = qJD(2) * qJD(3);
t97 = t89 * t99;
t58 = t97 + t101;
t73 = t89 * qJDD(2);
t98 = t86 * t99;
t59 = t73 - t98;
t27 = -t51 * qJD(4) + t88 * t58 + t85 * t59;
t78 = qJD(3) + qJD(4);
t48 = t78 * t51;
t116 = t27 - t48;
t102 = -g(3) + qJDD(1);
t82 = sin(pkin(7));
t83 = cos(pkin(7));
t64 = -t83 * g(1) - t82 * g(2);
t87 = sin(qJ(2));
t90 = cos(qJ(2));
t46 = t87 * t102 + t90 * t64;
t92 = qJD(2) ^ 2;
t39 = -t92 * pkin(2) + qJDD(2) * pkin(5) + t46;
t63 = -t82 * g(1) + t83 * g(2);
t30 = t86 * t39 - t89 * t63;
t114 = -t30 + (-t58 + t97) * pkin(6);
t49 = t51 ^ 2;
t50 = t53 ^ 2;
t76 = t78 ^ 2;
t31 = t89 * t39 + t86 * t63;
t67 = qJD(3) * pkin(3) - pkin(6) * t103;
t80 = t89 ^ 2;
t75 = t80 * t92;
t16 = -pkin(3) * t75 + t59 * pkin(6) - qJD(3) * t67 + t31;
t68 = t86 * t92 * t89;
t100 = qJDD(3) + t68;
t93 = t100 * pkin(3) + t114;
t7 = t85 * t16 - t88 * t93;
t107 = t88 * t16;
t8 = t85 * t93 + t107;
t2 = -t88 * t7 + t85 * t8;
t113 = t86 * t2;
t112 = t78 * t85;
t111 = t78 * t88;
t45 = t90 * t102 - t87 * t64;
t38 = -qJDD(2) * pkin(2) - t92 * pkin(5) - t45;
t25 = -t59 * pkin(3) - pkin(6) * t75 + t67 * t103 + t38;
t110 = t85 * t25;
t34 = t37 + t77;
t109 = t85 * t34;
t108 = t86 * t100;
t106 = t88 * t25;
t105 = t88 * t34;
t104 = t89 * (qJDD(3) - t68);
t3 = t85 * t7 + t88 * t8;
t12 = t86 * t30 + t89 * t31;
t96 = t85 * t58 - t88 * t59;
t94 = (-qJD(4) + t78) * t53 - t96;
t91 = qJD(3) ^ 2;
t79 = t86 ^ 2;
t74 = t79 * t92;
t62 = t74 + t75;
t61 = (t79 + t80) * qJDD(2);
t60 = t73 - 0.2e1 * t98;
t57 = 0.2e1 * t97 + t101;
t44 = -t50 + t76;
t43 = t49 - t76;
t42 = -t50 - t76;
t41 = -t104 - t86 * (-t74 - t91);
t40 = t89 * (-t75 - t91) - t108;
t36 = t50 - t49;
t32 = -t76 - t49;
t28 = -t49 - t50;
t26 = -t53 * qJD(4) - t96;
t24 = -t85 * t42 - t105;
t23 = t88 * t42 - t109;
t22 = t27 + t48;
t17 = (qJD(4) + t78) * t53 + t96;
t14 = t88 * t32 - t118;
t13 = t85 * t32 + t117;
t11 = -t86 * t23 + t89 * t24;
t10 = t85 * t22 + t88 * t94;
t9 = -t88 * t22 + t85 * t94;
t5 = -t86 * t13 + t89 * t14;
t4 = t89 * t10 - t86 * t9;
t1 = t89 * t3 - t113;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t102, 0, 0, 0, 0, 0, 0, t90 * qJDD(2) - t87 * t92, -t87 * qJDD(2) - t90 * t92, 0, t90 * t45 + t87 * t46, 0, 0, 0, 0, 0, 0, t87 * t40 + t90 * t60, t87 * t41 - t90 * t57, t87 * t61 + t90 * t62, t87 * t12 - t90 * t38, 0, 0, 0, 0, 0, 0, -t90 * t17 + t87 * t5, t87 * t11 - t116 * t90, -t90 * t28 + t87 * t4, t87 * t1 - t90 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t45, -t46, 0, 0, (t58 + t97) * t86, t89 * t57 + t86 * t60, t108 + t89 * (-t74 + t91), (t59 - t98) * t89, t86 * (t75 - t91) + t104, 0, pkin(2) * t60 + pkin(5) * t40 - t89 * t38, -pkin(2) * t57 + pkin(5) * t41 + t86 * t38, pkin(2) * t62 + pkin(5) * t61 + t12, -pkin(2) * t38 + pkin(5) * t12, t86 * (-t53 * t112 + t88 * t27) + t89 * (t53 * t111 + t85 * t27), t86 * (-t116 * t85 - t88 * t17) + t89 * (t116 * t88 - t85 * t17), t86 * (-t85 * t44 + t117) + t89 * (t88 * t44 + t118), t86 * (t51 * t111 - t85 * t26) + t89 * (t51 * t112 + t88 * t26), t86 * (t88 * t43 - t109) + t89 * (t85 * t43 + t105), (t86 * (-t51 * t88 + t53 * t85) + t89 * (-t51 * t85 - t53 * t88)) * t78, t86 * (-pkin(6) * t13 + t110) + t89 * (-pkin(3) * t17 + pkin(6) * t14 - t106) - pkin(2) * t17 + pkin(5) * t5, t86 * (-pkin(6) * t23 + t106) + t89 * (-pkin(3) * t116 + pkin(6) * t24 + t110) - pkin(2) * t116 + pkin(5) * t11, t86 * (-pkin(6) * t9 - t2) + t89 * (-pkin(3) * t28 + pkin(6) * t10 + t3) - pkin(2) * t28 + pkin(5) * t4, -pkin(6) * t113 + t89 * (-pkin(3) * t25 + pkin(6) * t3) - pkin(2) * t25 + pkin(5) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, t74 - t75, t101, t68, t73, qJDD(3), -t30, -t31, 0, 0, t37, t36, t22, -t37, t94, t77, pkin(3) * t13 - t7, -t107 - t85 * t114 + (-t100 * t85 + t23) * pkin(3), pkin(3) * t9, pkin(3) * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, t36, t22, -t37, t94, t77, -t7, -t8, 0, 0;];
tauJ_reg = t6;
