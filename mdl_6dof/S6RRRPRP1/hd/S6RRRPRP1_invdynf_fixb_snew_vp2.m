% Calculate vector of cutting forces with Newton-Euler
% S6RRRPRP1
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-05-07 07:28
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRPRP1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP1_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP1_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP1_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP1_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP1_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP1_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP1_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 07:24:28
% EndTime: 2019-05-07 07:24:38
% DurationCPUTime: 3.89s
% Computational Cost: add. (49997->205), mult. (112256->262), div. (0->0), fcn. (82135->10), ass. (0->96)
t100 = cos(pkin(10));
t102 = sin(qJ(3));
t106 = cos(qJ(3));
t103 = sin(qJ(2));
t107 = cos(qJ(2));
t128 = qJD(1) * qJD(2);
t109 = qJD(1) ^ 2;
t104 = sin(qJ(1));
t108 = cos(qJ(1));
t118 = -g(1) * t108 - g(2) * t104;
t86 = -pkin(1) * t109 + qJDD(1) * pkin(7) + t118;
t132 = t103 * t86;
t137 = pkin(2) * t109;
t89 = qJDD(1) * t103 + t107 * t128;
t58 = qJDD(2) * pkin(2) - pkin(8) * t89 - t132 + (pkin(8) * t128 + t103 * t137 - g(3)) * t107;
t122 = -g(3) * t103 + t107 * t86;
t90 = qJDD(1) * t107 - t103 * t128;
t130 = qJD(1) * t103;
t93 = qJD(2) * pkin(2) - pkin(8) * t130;
t98 = t107 ^ 2;
t59 = pkin(8) * t90 - qJD(2) * t93 - t98 * t137 + t122;
t119 = -t102 * t59 + t106 * t58;
t83 = (-t102 * t103 + t106 * t107) * qJD(1);
t66 = qJD(3) * t83 + t102 * t90 + t106 * t89;
t84 = (t102 * t107 + t103 * t106) * qJD(1);
t96 = qJDD(2) + qJDD(3);
t97 = qJD(2) + qJD(3);
t26 = (t83 * t97 - t66) * qJ(4) + (t83 * t84 + t96) * pkin(3) + t119;
t133 = t102 * t58 + t106 * t59;
t65 = -qJD(3) * t84 - t102 * t89 + t106 * t90;
t79 = pkin(3) * t97 - qJ(4) * t84;
t82 = t83 ^ 2;
t29 = -pkin(3) * t82 + t65 * qJ(4) - t79 * t97 + t133;
t99 = sin(pkin(10));
t76 = t100 * t84 + t83 * t99;
t141 = -0.2e1 * qJD(4) * t76 + t100 * t26 - t99 * t29;
t123 = t104 * g(1) - t108 * g(2);
t116 = -qJDD(1) * pkin(1) - t123;
t113 = -pkin(2) * t90 + t93 * t130 + (-pkin(8) * t98 - pkin(7)) * t109 + t116;
t101 = sin(qJ(5));
t105 = cos(qJ(5));
t75 = t100 * t83 - t84 * t99;
t124 = 0.2e1 * qJD(4) * t75 + t100 * t29 + t99 * t26;
t54 = -pkin(4) * t75 - pkin(9) * t76;
t95 = t97 ^ 2;
t21 = -pkin(4) * t95 + pkin(9) * t96 + t54 * t75 + t124;
t111 = -t65 * pkin(3) - qJ(4) * t82 + t84 * t79 + qJDD(4) + t113;
t43 = t100 * t65 - t66 * t99;
t44 = t100 * t66 + t65 * t99;
t24 = (-t75 * t97 - t44) * pkin(9) + (t76 * t97 - t43) * pkin(4) + t111;
t120 = -t101 * t21 + t105 * t24;
t62 = -t101 * t76 + t105 * t97;
t33 = t62 * qJD(5) + t101 * t96 + t105 * t44;
t42 = qJDD(5) - t43;
t72 = qJD(5) - t75;
t47 = -mrSges(7,2) * t72 + t62 * mrSges(7,3);
t63 = t101 * t97 + t105 * t76;
t127 = m(7) * (-0.2e1 * qJD(6) * t63 + (t62 * t72 - t33) * qJ(6) + (t62 * t63 + t42) * pkin(5) + t120) + t72 * t47 + t42 * mrSges(7,1);
t45 = -mrSges(7,1) * t62 + mrSges(7,2) * t63;
t46 = -mrSges(6,1) * t62 + mrSges(6,2) * t63;
t48 = -mrSges(6,2) * t72 + t62 * mrSges(6,3);
t11 = m(6) * t120 + t42 * mrSges(6,1) + t72 * t48 + (-t46 - t45) * t63 + (-mrSges(6,3) - mrSges(7,3)) * t33 + t127;
t135 = t101 * t24 + t105 * t21;
t32 = -t63 * qJD(5) - t101 * t44 + t105 * t96;
t49 = pkin(5) * t72 - t63 * qJ(6);
t60 = t62 ^ 2;
t126 = m(7) * (-t60 * pkin(5) + t32 * qJ(6) + 0.2e1 * qJD(6) * t62 - t49 * t72 + t135) + t32 * mrSges(7,3) + t62 * t45;
t50 = mrSges(7,1) * t72 - t63 * mrSges(7,3);
t51 = mrSges(6,1) * t72 - t63 * mrSges(6,3);
t13 = m(6) * t135 + t32 * mrSges(6,3) + t62 * t46 + (-t51 - t50) * t72 + (-mrSges(6,2) - mrSges(7,2)) * t42 + t126;
t68 = -mrSges(5,2) * t97 + mrSges(5,3) * t75;
t69 = mrSges(5,1) * t97 - mrSges(5,3) * t76;
t115 = -m(5) * t111 + t43 * mrSges(5,1) - t44 * mrSges(5,2) - t101 * t13 - t105 * t11 + t75 * t68 - t76 * t69;
t78 = -mrSges(4,2) * t97 + mrSges(4,3) * t83;
t80 = mrSges(4,1) * t97 - mrSges(4,3) * t84;
t112 = -m(4) * t113 + t65 * mrSges(4,1) - t66 * mrSges(4,2) + t83 * t78 - t84 * t80 + t115;
t91 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t130;
t129 = qJD(1) * t107;
t92 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t129;
t140 = (t103 * t91 - t107 * t92) * qJD(1) + m(3) * (-pkin(7) * t109 + t116) - t90 * mrSges(3,1) + t89 * mrSges(3,2) - t112;
t20 = -pkin(4) * t96 - pkin(9) * t95 + t76 * t54 - t141;
t125 = m(7) * (-t32 * pkin(5) - t60 * qJ(6) + t63 * t49 + qJDD(6) + t20) + t33 * mrSges(7,2) + t63 * t50;
t139 = -m(6) * t20 - t33 * mrSges(6,2) + (t48 + t47) * t62 + (mrSges(6,1) + mrSges(7,1)) * t32 - t63 * t51 - t125;
t53 = -mrSges(5,1) * t75 + mrSges(5,2) * t76;
t14 = m(5) * t141 + t96 * mrSges(5,1) - t44 * mrSges(5,3) - t76 * t53 + t97 * t68 + t139;
t77 = -mrSges(4,1) * t83 + mrSges(4,2) * t84;
t9 = m(5) * t124 - t96 * mrSges(5,2) + t43 * mrSges(5,3) - t101 * t11 + t105 * t13 + t75 * t53 - t97 * t69;
t6 = m(4) * t119 + t96 * mrSges(4,1) - t66 * mrSges(4,3) + t100 * t14 - t84 * t77 + t97 * t78 + t99 * t9;
t7 = m(4) * t133 - t96 * mrSges(4,2) + t65 * mrSges(4,3) + t100 * t9 - t99 * t14 + t83 * t77 - t97 * t80;
t88 = (-mrSges(3,1) * t107 + mrSges(3,2) * t103) * qJD(1);
t4 = m(3) * (-g(3) * t107 - t132) - t89 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t88 * t130 + qJD(2) * t92 + t102 * t7 + t106 * t6;
t5 = m(3) * t122 - qJDD(2) * mrSges(3,2) + t90 * mrSges(3,3) - qJD(2) * t91 - t102 * t6 + t106 * t7 + t88 * t129;
t138 = t103 * t5 + t107 * t4;
t8 = m(2) * t123 + qJDD(1) * mrSges(2,1) - t109 * mrSges(2,2) - t140;
t1 = m(2) * t118 - t109 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t103 * t4 + t107 * t5;
t2 = [-m(1) * g(1) + t1 * t108 - t104 * t8, t1, t5, t7, t9, t13, -t42 * mrSges(7,2) - t72 * t50 + t126; -m(1) * g(2) + t1 * t104 + t108 * t8, t8, t4, t6, t14, t11, -t33 * mrSges(7,3) - t63 * t45 + t127; (-m(1) - m(2)) * g(3) + t138, -m(2) * g(3) + t138, t140, -t112, -t115, -t139, -t32 * mrSges(7,1) - t62 * t47 + t125;];
f_new  = t2;
