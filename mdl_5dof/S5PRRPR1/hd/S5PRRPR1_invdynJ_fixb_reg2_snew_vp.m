% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRRPR1
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRRPR1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR1_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR1_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR1_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:15:51
% EndTime: 2019-12-05 16:15:55
% DurationCPUTime: 0.93s
% Computational Cost: add. (3256->137), mult. (4831->209), div. (0->0), fcn. (3378->10), ass. (0->107)
t106 = (qJD(2) + qJD(3));
t103 = t106 ^ 2;
t104 = qJDD(2) + qJDD(3);
t113 = sin(qJ(3));
t116 = cos(qJ(3));
t114 = sin(qJ(2));
t117 = cos(qJ(2));
t109 = sin(pkin(8));
t111 = cos(pkin(8));
t88 = g(1) * t109 - g(2) * t111;
t89 = -g(1) * t111 - g(2) * t109;
t128 = -t114 * t89 + t117 * t88;
t123 = qJDD(2) * pkin(2) + t128;
t127 = -t114 * t88 - t117 * t89;
t66 = -qJD(2) ^ 2 * pkin(2) - t127;
t50 = t113 * t123 + t116 * t66;
t160 = -pkin(3) * t103 + qJ(4) * t104 + (2 * qJD(4) * t106) + t50;
t112 = sin(qJ(5));
t108 = sin(pkin(9));
t110 = cos(pkin(9));
t115 = cos(qJ(5));
t157 = -t108 * t112 + t110 * t115;
t75 = t157 * t106;
t125 = t108 * t115 + t110 * t112;
t77 = t125 * t106;
t64 = t77 * t75;
t153 = qJDD(5) + t64;
t159 = t112 * t153;
t158 = t115 * t153;
t107 = -g(3) + qJDD(1);
t100 = t110 * t107;
t35 = t108 * t160 - t100;
t36 = t108 * t107 + t110 * t160;
t15 = t108 * t35 + t110 * t36;
t156 = t103 * t110;
t105 = t110 ^ 2;
t119 = t108 ^ 2;
t155 = t105 + t119;
t154 = t100 + (pkin(4) * t156 - pkin(7) * t104 - t160) * t108;
t85 = t155 * t103;
t73 = t75 ^ 2;
t74 = t77 ^ 2;
t142 = t110 * t104;
t145 = t105 * t103;
t29 = -pkin(4) * t145 + pkin(7) * t142 + t36;
t10 = t112 * t29 - t115 * t154;
t11 = t112 * t154 + t115 * t29;
t6 = -t10 * t115 + t11 * t112;
t152 = t108 * t6;
t49 = -t113 * t66 + t116 * t123;
t45 = -pkin(3) * t104 - qJ(4) * t103 + qJDD(4) - t49;
t151 = -pkin(3) * t45 + qJ(4) * t15;
t32 = -pkin(4) * t142 + t45 + (-t103 * t119 - t145) * pkin(7);
t150 = t112 * t32;
t58 = qJDD(5) - t64;
t149 = t112 * t58;
t148 = t115 * t32;
t147 = t115 * t58;
t144 = t108 * t104;
t140 = t116 * t104;
t139 = t75 * qJD(5);
t138 = t77 * qJD(5);
t136 = qJD(5) * t112;
t135 = qJD(5) * t115;
t81 = t155 * t156;
t134 = pkin(3) * t142 - qJ(4) * t81 - t110 * t45;
t51 = t157 * t104;
t72 = t125 * t104;
t39 = t112 * t51 - t115 * t72;
t40 = t112 * t72 + t115 * t51;
t20 = -t108 * t39 + t110 * t40;
t53 = -t73 - t74;
t7 = t10 * t112 + t11 * t115;
t132 = t108 * (-pkin(7) * t39 - t6) + t110 * (-pkin(4) * t53 + pkin(7) * t40 + t7) - pkin(3) * t53 + qJ(4) * t20;
t118 = qJD(5) ^ 2;
t56 = -t118 - t73;
t37 = t112 * t56 + t158;
t38 = t115 * t56 - t159;
t18 = -t108 * t37 + t110 * t38;
t60 = -t51 + 0.2e1 * t138;
t131 = t108 * (-pkin(7) * t37 + t150) + t110 * (-pkin(4) * t60 + pkin(7) * t38 - t148) - pkin(3) * t60 + qJ(4) * t18;
t69 = -t74 - t118;
t47 = t115 * t69 - t149;
t48 = -t112 * t69 - t147;
t25 = -t108 * t47 + t110 * t48;
t62 = t72 + 0.2e1 * t139;
t130 = t108 * (-pkin(7) * t47 + t148) + t110 * (-pkin(4) * t62 + pkin(7) * t48 + t150) - pkin(3) * t62 + qJ(4) * t25;
t97 = t119 * t104;
t98 = t105 * t104;
t84 = t98 + t97;
t129 = pkin(3) * t85 + qJ(4) * t84 + t15;
t80 = t108 * t85;
t126 = -pkin(3) * t144 + qJ(4) * t80 + t108 * t45;
t2 = t110 * t7 - t152;
t124 = -pkin(7) * t152 + qJ(4) * t2 - pkin(3) * t32 + t110 * (-pkin(4) * t32 + pkin(7) * t7);
t86 = 0.2e1 * t108 * t142;
t68 = -t74 + t118;
t67 = t73 - t118;
t63 = t72 + t139;
t61 = t51 - t138;
t30 = (t108 * (t112 * t77 + t115 * t75) + t110 * (t112 * t75 - t115 * t77)) * qJD(5);
t27 = t108 * (t115 * t63 - t136 * t77) + t110 * (t112 * t63 + t135 * t77);
t26 = t108 * (-t112 * t61 - t135 * t75) + t110 * (t115 * t61 - t136 * t75);
t24 = t108 * (-t112 * t68 + t158) + t110 * (t115 * t68 + t159);
t23 = t108 * (t115 * t67 - t149) + t110 * (t112 * t67 + t147);
t19 = t108 * (-t112 * t62 - t115 * t60) + t110 * (-t112 * t60 + t115 * t62);
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t107, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107, 0, 0, 0, 0, 0, 0, 0, 0, 0, t108 * t36 - t110 * t35, 0, 0, 0, 0, 0, 0, t108 * t38 + t110 * t37, t108 * t48 + t110 * t47, t108 * t40 + t110 * t39, t108 * t7 + t110 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t128, t127, 0, 0, 0, 0, 0, 0, 0, t104, pkin(2) * (-t103 * t113 + t140) + t49, pkin(2) * (-t103 * t116 - t104 * t113) - t50, 0, pkin(2) * (t113 * t50 + t116 * t49), t97, t86, 0, t98, 0, 0, pkin(2) * (t110 * t140 - t113 * t81) + t134, pkin(2) * (-t108 * t140 + t113 * t80) + t126, pkin(2) * (t113 * t84 + t116 * t85) + t129, pkin(2) * (t113 * t15 - t116 * t45) + t151, t27, t19, t24, t26, t23, t30, pkin(2) * (t113 * t18 - t116 * t60) + t131, pkin(2) * (t113 * t25 - t116 * t62) + t130, pkin(2) * (t113 * t20 - t116 * t53) + t132, pkin(2) * (t113 * t2 - t116 * t32) + t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104, t49, -t50, 0, 0, t97, t86, 0, t98, 0, 0, t134, t126, t129, t151, t27, t19, t24, t26, t23, t30, t131, t130, t132, t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t142, t144, -t85, t45, 0, 0, 0, 0, 0, 0, t60, t62, t53, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, t74 - t73, t72, t64, t51, qJDD(5), -t10, -t11, 0, 0;];
tauJ_reg = t1;
