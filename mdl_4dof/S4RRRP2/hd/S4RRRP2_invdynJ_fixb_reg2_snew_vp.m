% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
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
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4RRRP2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP2_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_invdynJ_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:13:09
% EndTime: 2019-12-31 17:13:11
% DurationCPUTime: 0.40s
% Computational Cost: add. (1172->105), mult. (1678->125), div. (0->0), fcn. (842->6), ass. (0->81)
t76 = qJD(1) + qJD(2);
t102 = (qJD(4) * t76);
t118 = 2 * t102;
t74 = t76 ^ 2;
t80 = sin(qJ(3));
t83 = cos(qJ(3));
t62 = t83 * t74 * t80;
t56 = qJDD(3) + t62;
t117 = pkin(3) * t56;
t82 = sin(qJ(1));
t85 = cos(qJ(1));
t96 = t82 * g(1) - t85 * g(2);
t53 = qJDD(1) * pkin(1) + t96;
t92 = t85 * g(1) + t82 * g(2);
t54 = -qJD(1) ^ 2 * pkin(1) - t92;
t81 = sin(qJ(2));
t84 = cos(qJ(2));
t30 = t81 * t53 + t84 * t54;
t75 = qJDD(1) + qJDD(2);
t25 = -t74 * pkin(2) + t75 * pkin(6) + t30;
t112 = t80 * t25;
t116 = t83 * g(3);
t18 = t112 + t116;
t19 = -t80 * g(3) + t83 * t25;
t8 = t80 * t18 + t83 * t19;
t69 = t83 * t75;
t103 = qJD(3) * t76;
t98 = t80 * t103;
t46 = t69 - t98;
t104 = qJ(4) * t80;
t55 = qJD(3) * pkin(3) - t76 * t104;
t88 = t46 * qJ(4) - qJD(3) * t55 + t83 * t118 + t19;
t29 = t84 * t53 - t81 * t54;
t24 = -t75 * pkin(2) - t74 * pkin(6) - t29;
t115 = -pkin(2) * t24 + pkin(6) * t8;
t77 = t80 ^ 2;
t114 = t77 * t74;
t78 = t83 ^ 2;
t113 = t78 * t74;
t111 = t80 * t56;
t68 = t80 * t75;
t57 = qJDD(3) - t62;
t110 = t83 * t57;
t86 = qJD(3) ^ 2;
t60 = -t86 - t113;
t37 = t83 * t60 - t111;
t47 = t69 - 0.2e1 * t98;
t109 = pkin(2) * t47 + pkin(6) * t37;
t59 = -t86 - t114;
t38 = -t80 * t59 - t110;
t97 = t83 * t103;
t44 = t68 + 0.2e1 * t97;
t108 = -pkin(2) * t44 + pkin(6) * t38;
t106 = t77 + t78;
t50 = t106 * t75;
t51 = t106 * t74;
t107 = pkin(2) * t51 + pkin(6) * t50;
t105 = qJ(4) * t75;
t101 = t80 * t24 + t108;
t100 = -t83 * t24 + t109;
t45 = t68 + t97;
t87 = t116 + (t45 - t97) * qJ(4) - t117;
t95 = t80 * ((t118 + t25 + t105) * t80 + t87) + t83 * (t83 * t105 + (t51 - t113) * pkin(3) + t88) + t107;
t13 = t80 * t76 * t55 - t46 * pkin(3) - qJ(4) * t113 + qJDD(4) + t24;
t94 = t80 * (-qJ(4) * t59 + t13) + t83 * (-pkin(3) * t44 - qJ(4) * t57) + t108;
t93 = t107 + t8;
t10 = -0.2e1 * t80 * t102 - t112 - t87;
t11 = -pkin(3) * t113 + t88;
t2 = -t80 * t10 + t83 * t11;
t90 = pkin(6) * t2 - t10 * t104 - pkin(2) * t13 + t83 * (-pkin(3) * t13 + qJ(4) * t11);
t89 = -t56 * t104 + t109 + t83 * (pkin(3) * t47 + qJ(4) * t60 - t13);
t52 = (t77 - t78) * t74;
t36 = t111 + t83 * (t86 - t114);
t35 = t80 * (-t86 + t113) + t110;
t32 = (t45 + t97) * t80;
t31 = (t46 - t98) * t83;
t27 = pkin(1) * (t81 * t50 + t84 * t51);
t26 = t83 * t44 + t80 * t47;
t17 = pkin(1) * (t81 * t38 - t84 * t44);
t16 = pkin(1) * (t81 * t37 + t84 * t47);
t1 = [0, 0, 0, 0, 0, qJDD(1), t96, t92, 0, 0, 0, 0, 0, 0, 0, t75, pkin(1) * (-t81 * t74 + t84 * t75) + t29, pkin(1) * (-t84 * t74 - t81 * t75) - t30, 0, pkin(1) * (t84 * t29 + t81 * t30), t32, t26, t36, t31, t35, 0, t16 + t100, t17 + t101, t27 + t93, pkin(1) * (-t84 * t24 + t81 * t8) + t115, t32, t26, t36, t31, t35, 0, t16 + t89, t17 + t94, t27 + t95, pkin(1) * (-t84 * t13 + t81 * t2) + t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t29, -t30, 0, 0, t32, t26, t36, t31, t35, 0, t100, t101, t93, t115, t32, t26, t36, t31, t35, 0, t89, t94, t95, t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, t52, t68, t62, t69, qJDD(3), -t18, -t19, 0, 0, -t62, t52, t68, t62, t69, qJDD(3), t10 + t117, (t59 + t113) * pkin(3) - t88, -pkin(3) * t68, pkin(3) * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, t44, -t51, t13;];
tauJ_reg = t1;
