% Calculate vector of cutting forces with Newton-Euler
% S6RRPRRP10
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% f_new [3x7]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 18:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRRP10_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP10_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP10_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP10_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP10_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP10_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP10_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP10_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP10_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 18:39:18
% EndTime: 2019-05-06 18:39:33
% DurationCPUTime: 6.29s
% Computational Cost: add. (81434->211), mult. (184713->281), div. (0->0), fcn. (147739->12), ass. (0->106)
t108 = sin(qJ(4));
t111 = cos(qJ(4));
t103 = sin(pkin(11));
t105 = cos(pkin(11));
t104 = sin(pkin(6));
t109 = sin(qJ(2));
t112 = cos(qJ(2));
t132 = qJD(1) * t112;
t106 = cos(pkin(6));
t114 = qJD(1) ^ 2;
t110 = sin(qJ(1));
t113 = cos(qJ(1));
t124 = t110 * g(1) - t113 * g(2);
t88 = t114 * t104 * pkin(8) + qJDD(1) * pkin(1) + t124;
t136 = t106 * t88;
t121 = -t113 * g(1) - t110 * g(2);
t130 = qJDD(1) * t104;
t89 = -t114 * pkin(1) + pkin(8) * t130 + t121;
t137 = t109 * t136 + t112 * t89;
t133 = qJD(1) * t104;
t90 = (-pkin(2) * t112 - qJ(3) * t109) * t133;
t100 = t106 * qJD(1) + qJD(2);
t98 = t100 ^ 2;
t99 = t106 * qJDD(1) + qJDD(2);
t58 = -t98 * pkin(2) + t99 * qJ(3) + (-g(3) * t109 + t90 * t132) * t104 + t137;
t144 = t106 * g(3);
t92 = (qJD(2) * t132 + qJDD(1) * t109) * t104;
t126 = t109 * t133;
t93 = -qJD(2) * t126 + t112 * t130;
t59 = -t93 * pkin(2) - t144 - t92 * qJ(3) + (-t88 + (pkin(2) * t109 - qJ(3) * t112) * t100 * qJD(1)) * t104;
t82 = t103 * t100 + t105 * t126;
t122 = -0.2e1 * qJD(3) * t82 - t103 * t58 + t105 * t59;
t125 = t104 * t132;
t73 = t103 * t99 + t105 * t92;
t81 = t105 * t100 - t103 * t126;
t27 = (-t81 * t125 - t73) * pkin(9) + (t81 * t82 - t93) * pkin(3) + t122;
t127 = 0.2e1 * qJD(3) * t81 + t103 * t59 + t105 * t58;
t72 = -t103 * t92 + t105 * t99;
t74 = -pkin(3) * t125 - t82 * pkin(9);
t79 = t81 ^ 2;
t30 = -t79 * pkin(3) + t72 * pkin(9) + t74 * t125 + t127;
t123 = -t108 * t30 + t111 * t27;
t67 = -t108 * t82 + t111 * t81;
t68 = t108 * t81 + t111 * t82;
t53 = -t67 * pkin(4) - t68 * pkin(10);
t85 = qJDD(4) - t93;
t96 = qJD(4) - t125;
t95 = t96 ^ 2;
t22 = -t85 * pkin(4) - t95 * pkin(10) + t68 * t53 - t123;
t107 = sin(qJ(5));
t145 = cos(qJ(5));
t46 = t67 * qJD(4) + t108 * t72 + t111 * t73;
t61 = t107 * t96 + t145 * t68;
t32 = t61 * qJD(5) + t107 * t46 - t145 * t85;
t60 = t107 * t68 - t145 * t96;
t33 = -t60 * qJD(5) + t107 * t85 + t145 * t46;
t66 = qJD(5) - t67;
t47 = -t60 * mrSges(7,2) + t66 * mrSges(7,3);
t128 = m(7) * (-0.2e1 * qJD(6) * t61 + (t60 * t66 - t33) * qJ(6) + (t61 * t66 + t32) * pkin(5) + t22) + t32 * mrSges(7,1) + t60 * t47;
t48 = -t66 * mrSges(6,2) - t60 * mrSges(6,3);
t49 = t66 * mrSges(6,1) - t61 * mrSges(6,3);
t50 = -t66 * mrSges(7,1) + t61 * mrSges(7,2);
t147 = m(6) * t22 + t32 * mrSges(6,1) + (t49 - t50) * t61 + (mrSges(6,2) - mrSges(7,3)) * t33 + t60 * t48 + t128;
t140 = t108 * t27 + t111 * t30;
t23 = -t95 * pkin(4) + t85 * pkin(10) + t67 * t53 + t140;
t134 = t104 * t112;
t120 = -g(3) * t134 - t109 * t89 + t112 * t136;
t57 = -t99 * pkin(2) - t98 * qJ(3) + t90 * t126 + qJDD(3) - t120;
t116 = -t72 * pkin(3) - t79 * pkin(9) + t82 * t74 + t57;
t45 = -t68 * qJD(4) - t108 * t73 + t111 * t72;
t25 = (-t67 * t96 - t46) * pkin(10) + (t68 * t96 - t45) * pkin(4) + t116;
t119 = -t107 * t23 + t145 * t25;
t39 = t60 * pkin(5) - t61 * qJ(6);
t44 = qJDD(5) - t45;
t65 = t66 ^ 2;
t146 = m(7) * (-t44 * pkin(5) - t65 * qJ(6) + t61 * t39 + qJDD(6) - t119);
t142 = -mrSges(6,3) - mrSges(7,2);
t141 = t107 * t25 + t145 * t23;
t40 = t60 * mrSges(7,1) - t61 * mrSges(7,3);
t139 = -t60 * mrSges(6,1) - t61 * mrSges(6,2) - t40;
t135 = t104 * t109;
t129 = m(7) * (-t65 * pkin(5) + t44 * qJ(6) + 0.2e1 * qJD(6) * t66 - t60 * t39 + t141) + t66 * t50 + t44 * mrSges(7,3);
t14 = m(6) * t141 - t44 * mrSges(6,2) + t139 * t60 + t142 * t32 - t66 * t49 + t129;
t16 = m(6) * t119 - t146 + (t48 + t47) * t66 + t139 * t61 + (mrSges(6,1) + mrSges(7,1)) * t44 + t142 * t33;
t62 = -t96 * mrSges(5,2) + t67 * mrSges(5,3);
t63 = t96 * mrSges(5,1) - t68 * mrSges(5,3);
t118 = -m(5) * t116 + t45 * mrSges(5,1) - t46 * mrSges(5,2) - t107 * t14 - t145 * t16 + t67 * t62 - t68 * t63;
t70 = mrSges(4,2) * t125 + t81 * mrSges(4,3);
t71 = -mrSges(4,1) * t125 - t82 * mrSges(4,3);
t115 = m(4) * t57 - t72 * mrSges(4,1) + t73 * mrSges(4,2) - t81 * t70 + t82 * t71 - t118;
t87 = -t100 * mrSges(3,2) + mrSges(3,3) * t125;
t91 = (-mrSges(3,1) * t112 + mrSges(3,2) * t109) * t133;
t10 = m(3) * t120 + t99 * mrSges(3,1) - t92 * mrSges(3,3) + t100 * t87 - t91 * t126 - t115;
t52 = -t67 * mrSges(5,1) + t68 * mrSges(5,2);
t11 = m(5) * t140 - t85 * mrSges(5,2) + t45 * mrSges(5,3) - t107 * t16 + t145 * t14 + t67 * t52 - t96 * t63;
t12 = m(5) * t123 + t85 * mrSges(5,1) - t46 * mrSges(5,3) - t68 * t52 + t96 * t62 - t147;
t69 = -t81 * mrSges(4,1) + t82 * mrSges(4,2);
t7 = m(4) * t122 - t93 * mrSges(4,1) - t73 * mrSges(4,3) + t108 * t11 + t111 * t12 - t70 * t125 - t82 * t69;
t8 = m(4) * t127 + t93 * mrSges(4,2) + t72 * mrSges(4,3) - t108 * t12 + t111 * t11 + t71 * t125 + t81 * t69;
t86 = t100 * mrSges(3,1) - mrSges(3,3) * t126;
t4 = m(3) * (-g(3) * t135 + t137) + t93 * mrSges(3,3) - t99 * mrSges(3,2) + t91 * t125 - t100 * t86 + t105 * t8 - t103 * t7;
t6 = m(3) * (-t104 * t88 - t144) + t92 * mrSges(3,2) - t93 * mrSges(3,1) + t103 * t8 + t105 * t7 + (t109 * t86 - t112 * t87) * t133;
t131 = t10 * t134 + t106 * t6 + t4 * t135;
t2 = m(2) * t121 - t114 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t109 * t10 + t112 * t4;
t1 = m(2) * t124 + qJDD(1) * mrSges(2,1) - t114 * mrSges(2,2) - t104 * t6 + (t112 * t10 + t109 * t4) * t106;
t3 = [-m(1) * g(1) - t110 * t1 + t113 * t2, t2, t4, t8, t11, t14, -t32 * mrSges(7,2) - t60 * t40 + t129; -m(1) * g(2) + t113 * t1 + t110 * t2, t1, t10, t7, t12, t16, -t33 * mrSges(7,3) - t61 * t50 + t128; (-m(1) - m(2)) * g(3) + t131, -m(2) * g(3) + t131, t6, t115, -t118, t147, -t44 * mrSges(7,1) + t33 * mrSges(7,2) + t61 * t40 - t66 * t47 + t146;];
f_new  = t3;
