% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRRPR2
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
% Datum: 2019-12-05 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRRPR2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR2_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:17:35
% EndTime: 2019-12-05 16:17:38
% DurationCPUTime: 0.76s
% Computational Cost: add. (2770->126), mult. (4052->182), div. (0->0), fcn. (2723->10), ass. (0->106)
t103 = sin(pkin(9));
t107 = sin(qJ(5));
t105 = cos(pkin(9));
t99 = qJD(2) + qJD(3);
t138 = t105 * t99;
t84 = -qJD(5) + t138;
t127 = t107 * t84 * t99;
t110 = cos(qJ(5));
t96 = qJDD(2) + qJDD(3);
t134 = t110 * t96;
t131 = qJD(5) * t99;
t78 = t107 * t103 * t131;
t115 = t103 * (t127 + t134) - t78;
t149 = t103 * t115;
t95 = t99 ^ 2;
t97 = t103 ^ 2;
t98 = t105 ^ 2;
t73 = (t97 + t98) * t95;
t145 = 2 * qJD(4);
t108 = sin(qJ(3));
t111 = cos(qJ(3));
t109 = sin(qJ(2));
t112 = cos(qJ(2));
t104 = sin(pkin(8));
t106 = cos(pkin(8));
t80 = t104 * g(1) - t106 * g(2);
t82 = -t106 * g(1) - t104 * g(2);
t123 = -t109 * t82 + t112 * t80;
t56 = qJDD(2) * pkin(2) + t123;
t120 = -t109 * t80 - t112 * t82;
t57 = -qJD(2) ^ 2 * pkin(2) - t120;
t44 = t108 * t56 + t111 * t57;
t40 = -t95 * pkin(3) + t96 * qJ(4) + t44;
t148 = t99 * t145 + t40;
t121 = -t105 * pkin(4) - t103 * pkin(7);
t102 = -g(3) + qJDD(1);
t92 = t105 * t102;
t29 = t148 * t103 - t92;
t30 = t103 * t102 + t148 * t105;
t12 = t103 * t29 + t105 * t30;
t136 = t107 * t96;
t146 = t103 * (t110 * t131 + t136);
t79 = t84 ^ 2;
t70 = t121 * t99;
t15 = -t92 + (t40 + (t145 + t70) * t99) * t103;
t16 = t70 * t138 + t30;
t43 = -t108 * t57 + t111 * t56;
t39 = -t96 * pkin(3) - t95 * qJ(4) + qJDD(4) - t43;
t31 = t121 * t96 + t39;
t10 = t107 * t31 + t110 * t16;
t9 = t107 * t16 - t110 * t31;
t6 = t110 * t10 + t107 * t9;
t2 = t103 * t15 + t105 * t6;
t5 = t107 * t10 - t110 * t9;
t144 = -pkin(3) * t5 + qJ(4) * t2;
t142 = t97 * t95;
t141 = -pkin(3) * t39 + qJ(4) * t12;
t140 = t103 * t96;
t139 = t105 * t96;
t126 = t107 * t142;
t76 = t110 * t126;
t81 = -qJDD(5) + t139;
t59 = -t76 + t81;
t137 = t107 * t59;
t60 = -t76 - t81;
t135 = t110 * t60;
t133 = t110 * t99;
t132 = t111 * t96;
t68 = t105 * t73;
t129 = pkin(3) * t139 - qJ(4) * t68 - t105 * t39;
t128 = t110 ^ 2 * t142;
t58 = -t128 - t79;
t47 = -t107 * t58 + t110 * t59;
t21 = t105 * t47 + t149;
t46 = t110 * t58 + t137;
t125 = t103 * (-pkin(7) * t46 + t110 * t15) + t105 * (-pkin(4) * t46 + t10) - pkin(3) * t46 + qJ(4) * t21;
t83 = t107 ^ 2 * t142;
t61 = -t83 - t79;
t50 = -t107 * t60 + t110 * t61;
t66 = t103 * t84 * t133;
t52 = t66 - t146;
t28 = -t103 * t52 + t105 * t50;
t49 = t107 * t61 + t135;
t124 = t103 * (-pkin(7) * t49 + t107 * t15) + t105 * (-pkin(4) * t49 + t9) - pkin(3) * t49 + qJ(4) * t28;
t89 = t97 * t96;
t90 = t98 * t96;
t72 = t90 + t89;
t122 = pkin(3) * t73 + qJ(4) * t72 + t12;
t67 = t103 * t73;
t119 = -pkin(3) * t140 + qJ(4) * t67 + t103 * t39;
t51 = t66 + t146;
t53 = -t78 + (-t127 + t134) * t103;
t34 = t107 * t53 - t110 * t51;
t62 = t83 + t128;
t19 = -t103 * t62 + t105 * t34;
t33 = -t107 * t51 - t110 * t53;
t118 = qJ(4) * t19 - t103 * t5 + (-pkin(3) + t121) * t33;
t77 = 0.2e1 * t103 * t139;
t74 = t105 * t81;
t63 = -t83 + t128;
t42 = (-t105 * t126 + t149) * t110;
t41 = t105 * t76 + (t136 + (qJD(5) - t84) * t133) * t107 * t97;
t27 = t103 * (t110 * (t83 - t79) + t137) + t105 * t51;
t26 = t103 * (t135 - t107 * (t79 - t128)) - t105 * t53;
t18 = t103 * (-t107 * t115 + t110 * t52) - t105 * t63;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t102, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103 * t30 - t105 * t29, 0, 0, 0, 0, 0, 0, t103 * t50 + t105 * t52, t103 * t47 - t105 * t115, t103 * t34 + t105 * t62, t103 * t6 - t105 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t123, t120, 0, 0, 0, 0, 0, 0, 0, t96, pkin(2) * (-t108 * t95 + t132) + t43, pkin(2) * (-t108 * t96 - t111 * t95) - t44, 0, pkin(2) * (t108 * t44 + t111 * t43), t89, t77, 0, t90, 0, 0, pkin(2) * (t105 * t132 - t108 * t68) + t129, pkin(2) * (-t103 * t132 + t108 * t67) + t119, pkin(2) * (t108 * t72 + t111 * t73) + t122, pkin(2) * (t108 * t12 - t111 * t39) + t141, t42, t18, t26, t41, t27, t74, pkin(2) * (t108 * t28 - t111 * t49) + t124, pkin(2) * (t108 * t21 - t111 * t46) + t125, pkin(2) * (t108 * t19 - t111 * t33) + t118, pkin(2) * t108 * t2 + (-pkin(2) * t111 + t121) * t5 + t144; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, t43, -t44, 0, 0, t89, t77, 0, t90, 0, 0, t129, t119, t122, t141, t42, t18, t26, t41, t27, t74, t124, t125, t118, t121 * t5 + t144; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t139, t140, -t73, t39, 0, 0, 0, 0, 0, 0, t49, t46, t33, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, t63, t53, -t76, -t51, -t81, -t9, -t10, 0, 0;];
tauJ_reg = t1;
