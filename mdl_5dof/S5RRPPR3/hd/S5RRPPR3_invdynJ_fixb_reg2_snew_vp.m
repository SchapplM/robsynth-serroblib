% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRPPR3
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRPPR3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR3_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR3_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR3_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:26:36
% EndTime: 2019-12-31 19:26:37
% DurationCPUTime: 0.53s
% Computational Cost: add. (1820->104), mult. (2454->146), div. (0->0), fcn. (1350->8), ass. (0->84)
t79 = (qJD(1) + qJD(2));
t77 = t79 ^ 2;
t78 = qJDD(1) + qJDD(2);
t83 = sin(pkin(8));
t84 = cos(pkin(8));
t55 = t84 * t77 + t83 * t78;
t58 = t83 * t77 - t84 * t78;
t86 = sin(qJ(2));
t89 = cos(qJ(2));
t124 = pkin(1) * (t86 * t55 + t89 * t58);
t123 = pkin(1) * (t89 * t55 - t86 * t58);
t122 = pkin(3) + pkin(7);
t121 = pkin(2) * t55;
t120 = pkin(2) * t58;
t82 = -g(3) + qJDD(3);
t85 = sin(qJ(5));
t88 = cos(qJ(5));
t87 = sin(qJ(1));
t90 = cos(qJ(1));
t103 = t87 * g(1) - t90 * g(2);
t61 = qJDD(1) * pkin(1) + t103;
t99 = t90 * g(1) + t87 * g(2);
t62 = -qJD(1) ^ 2 * pkin(1) - t99;
t34 = t89 * t61 - t86 * t62;
t32 = t78 * pkin(2) + t34;
t35 = t86 * t61 + t89 * t62;
t33 = -t77 * pkin(2) + t35;
t102 = -t84 * t32 + t83 * t33;
t117 = t78 * pkin(3);
t97 = qJDD(4) + t102 - t117;
t19 = -t77 * qJ(4) + t97;
t92 = -t78 * pkin(7) + t19;
t10 = t85 * t82 - t88 * t92;
t11 = t88 * t82 + t85 * t92;
t4 = -t88 * t10 + t85 * t11;
t80 = t85 ^ 2;
t116 = t80 * t77;
t81 = t88 ^ 2;
t115 = t81 * t77;
t106 = t85 * t77 * t88;
t63 = qJDD(5) + t106;
t113 = t85 * t63;
t112 = t85 * t78;
t64 = qJDD(5) - t106;
t111 = t88 * t64;
t21 = t83 * t32 + t84 * t33;
t110 = t80 + t81;
t109 = t78 * qJ(4);
t108 = qJD(5) * t79;
t105 = (2 * qJD(4) * t79) + t21;
t17 = -t77 * pkin(3) + t105 + t109;
t6 = t83 * t17 - t84 * t19;
t107 = pkin(2) * t6 - pkin(3) * t19 + qJ(4) * t17;
t104 = -t21 - t121;
t101 = -t102 - t120;
t100 = t105 + 0.2e1 * t109 + t121;
t67 = t88 * t78;
t53 = -0.2e1 * t85 * t108 + t67;
t14 = -t77 * pkin(7) + t17;
t2 = t83 * t14 - t84 * t4;
t98 = pkin(2) * t2 + qJ(4) * t14 - t122 * t4;
t91 = qJD(5) ^ 2;
t66 = -t91 - t115;
t41 = t88 * t66 - t113;
t25 = -t84 * t41 + t83 * t53;
t96 = pkin(2) * t25 + qJ(4) * t53 - t122 * t41 + t88 * t14;
t52 = 0.2e1 * t88 * t108 + t112;
t65 = -t91 - t116;
t40 = t85 * t65 + t111;
t24 = -t84 * t40 + t83 * t52;
t95 = pkin(2) * t24 + qJ(4) * t52 - t122 * t40 + t85 * t14;
t94 = -t117 + t97 + t120;
t59 = t110 * t78;
t60 = t110 * t77;
t31 = t84 * t59 - t83 * t60;
t93 = pkin(2) * t31 - qJ(4) * t60 + t122 * t59 - t4;
t43 = t111 - t85 * (t91 - t115);
t42 = t88 * (-t91 + t116) - t113;
t37 = t53 * t88;
t36 = t52 * t85;
t29 = -t88 * t52 - t85 * t53;
t8 = -t102 * t84 + t83 * t21;
t7 = pkin(2) * t8;
t1 = [0, 0, 0, 0, 0, qJDD(1), t103, t99, 0, 0, 0, 0, 0, 0, 0, t78, pkin(1) * (-t86 * t77 + t89 * t78) + t34, pkin(1) * (-t89 * t77 - t86 * t78) - t35, 0, pkin(1) * (t89 * t34 + t86 * t35), 0, 0, 0, 0, 0, t78, t101 - t124, t104 - t123, 0, pkin(1) * (t86 * (t102 * t83 + t84 * t21) + t89 * t8) + t7, t78, 0, 0, 0, 0, 0, 0, t94 + t124, t100 + t123, pkin(1) * (t86 * (t84 * t17 + t83 * t19) + t89 * t6) + t107, t37, t29, t43, t36, t42, 0, pkin(1) * (t86 * (t83 * t40 + t84 * t52) + t89 * t24) + t95, pkin(1) * (t86 * (t83 * t41 + t84 * t53) + t89 * t25) + t96, pkin(1) * (t86 * (-t83 * t59 - t84 * t60) + t89 * t31) + t93, pkin(1) * (t86 * (t84 * t14 + t83 * t4) + t89 * t2) + t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, t34, -t35, 0, 0, 0, 0, 0, 0, 0, t78, t101, t104, 0, t7, t78, 0, 0, 0, 0, 0, 0, t94, t100, t107, t37, t29, t43, t36, t42, 0, t95, t96, t93, t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, 0, 0, 0, 0, 0, 0, -t85 * t64 + t88 * t65, -t88 * t63 - t85 * t66, 0, t85 * t10 + t88 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, -t77, t19, 0, 0, 0, 0, 0, 0, t40, t41, -t59, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, (-t80 + t81) * t77, t67, -t106, -t112, qJDD(5), -t10, -t11, 0, 0;];
tauJ_reg = t1;
