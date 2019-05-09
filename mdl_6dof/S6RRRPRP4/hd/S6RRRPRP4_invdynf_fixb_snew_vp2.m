% Calculate vector of cutting forces with Newton-Euler
% S6RRRPRP4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-05-07 07:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRPRP4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP4_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP4_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP4_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP4_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP4_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 07:44:18
% EndTime: 2019-05-07 07:44:24
% DurationCPUTime: 1.78s
% Computational Cost: add. (18137->199), mult. (37318->243), div. (0->0), fcn. (24932->8), ass. (0->93)
t147 = -2 * qJD(4);
t100 = cos(qJ(2));
t102 = qJD(1) ^ 2;
t101 = cos(qJ(1));
t99 = sin(qJ(1));
t120 = t99 * g(1) - t101 * g(2);
t115 = -qJDD(1) * pkin(1) - t120;
t98 = sin(qJ(2));
t125 = qJD(1) * t98;
t122 = qJD(1) * qJD(2);
t84 = t100 * qJDD(1) - t122 * t98;
t87 = qJD(2) * pkin(2) - pkin(8) * t125;
t95 = t100 ^ 2;
t107 = -t84 * pkin(2) + t87 * t125 + (-pkin(8) * t95 - pkin(7)) * t102 + t115;
t139 = cos(qJ(5));
t123 = qJD(1) * t100;
t140 = cos(qJ(3));
t97 = sin(qJ(3));
t76 = -t123 * t140 + t125 * t97;
t94 = qJD(2) + qJD(3);
t137 = t76 * t94;
t83 = t98 * qJDD(1) + t100 * t122;
t50 = -t76 * qJD(3) + t140 * t83 + t97 * t84;
t77 = (t100 * t97 + t140 * t98) * qJD(1);
t103 = (-t50 + t137) * qJ(4) + t107 + (t94 * pkin(3) + t147) * t77;
t49 = t77 * qJD(3) - t140 * t84 + t97 * t83;
t68 = t77 * pkin(4) - t94 * pkin(9);
t75 = t76 ^ 2;
t17 = -t75 * pkin(4) - t77 * t68 + (pkin(3) + pkin(9)) * t49 + t103;
t117 = -t101 * g(1) - t99 * g(2);
t79 = -t102 * pkin(1) + qJDD(1) * pkin(7) + t117;
t134 = t98 * t79;
t138 = pkin(2) * t102;
t39 = qJDD(2) * pkin(2) - t83 * pkin(8) - t134 + (pkin(8) * t122 + t138 * t98 - g(3)) * t100;
t119 = -t98 * g(3) + t100 * t79;
t40 = t84 * pkin(8) - qJD(2) * t87 - t138 * t95 + t119;
t118 = t140 * t39 - t97 * t40;
t58 = t76 * pkin(3) - t77 * qJ(4);
t92 = t94 ^ 2;
t93 = qJDD(2) + qJDD(3);
t25 = -t93 * pkin(3) - t92 * qJ(4) + t77 * t58 + qJDD(4) - t118;
t19 = (t76 * t77 - t93) * pkin(9) + (t50 + t137) * pkin(4) + t25;
t96 = sin(qJ(5));
t130 = t139 * t17 + t96 * t19;
t62 = -t139 * t76 + t96 * t94;
t63 = t139 * t94 + t96 * t76;
t36 = t62 * pkin(5) - t63 * qJ(6);
t47 = qJDD(5) + t50;
t74 = qJD(5) + t77;
t55 = -t74 * mrSges(7,1) + t63 * mrSges(7,2);
t73 = t74 ^ 2;
t121 = m(7) * (-t73 * pkin(5) + t47 * qJ(6) + 0.2e1 * qJD(6) * t74 - t62 * t36 + t130) + t74 * t55 + t47 * mrSges(7,3);
t37 = t62 * mrSges(7,1) - t63 * mrSges(7,3);
t128 = -t62 * mrSges(6,1) - t63 * mrSges(6,2) - t37;
t131 = -mrSges(6,3) - mrSges(7,2);
t28 = t63 * qJD(5) - t139 * t49 + t96 * t93;
t54 = t74 * mrSges(6,1) - t63 * mrSges(6,3);
t10 = m(6) * t130 - t47 * mrSges(6,2) + t128 * t62 + t131 * t28 - t74 * t54 + t121;
t112 = t139 * t19 - t96 * t17;
t141 = m(7) * (-t47 * pkin(5) - t73 * qJ(6) + t63 * t36 + qJDD(6) - t112);
t29 = -t62 * qJD(5) + t139 * t93 + t96 * t49;
t52 = -t62 * mrSges(7,2) + t74 * mrSges(7,3);
t53 = -t74 * mrSges(6,2) - t62 * mrSges(6,3);
t11 = m(6) * t112 - t141 + (t53 + t52) * t74 + t128 * t63 + (mrSges(6,1) + mrSges(7,1)) * t47 + t131 * t29;
t67 = t77 * mrSges(5,1) + t94 * mrSges(5,2);
t113 = -t139 * t10 + t96 * t11 - m(5) * (t49 * pkin(3) + t103) + t50 * mrSges(5,3) + t77 * t67;
t66 = t76 * mrSges(5,1) - t94 * mrSges(5,3);
t126 = -t94 * mrSges(4,2) - t76 * mrSges(4,3) - t66;
t133 = mrSges(4,1) - mrSges(5,2);
t65 = t94 * mrSges(4,1) - t77 * mrSges(4,3);
t105 = m(4) * t107 + t50 * mrSges(4,2) + t126 * t76 + t133 * t49 + t77 * t65 - t113;
t85 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t125;
t86 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t123;
t146 = t105 - (t100 * t86 - t98 * t85) * qJD(1) - t84 * mrSges(3,1) + t83 * mrSges(3,2) + m(3) * (-t102 * pkin(7) + t115);
t129 = t140 * t40 + t97 * t39;
t145 = t92 * pkin(3) - t93 * qJ(4) + t94 * t147 + t76 * t58 - t129;
t110 = m(5) * t25 + t96 * t10 + t11 * t139;
t60 = -t76 * mrSges(5,2) - t77 * mrSges(5,3);
t127 = -t76 * mrSges(4,1) - t77 * mrSges(4,2) - t60;
t132 = -mrSges(4,3) - mrSges(5,1);
t7 = m(4) * t118 + t126 * t94 + t127 * t77 + t132 * t50 + t133 * t93 - t110;
t108 = -t49 * pkin(4) - t75 * pkin(9) + t94 * t68 - t145;
t111 = -t29 * mrSges(7,3) - t63 * t55 + m(7) * (-0.2e1 * qJD(6) * t63 + (t62 * t74 - t29) * qJ(6) + (t63 * t74 + t28) * pkin(5) + t108) + t28 * mrSges(7,1) + t62 * t52;
t106 = m(6) * t108 + t28 * mrSges(6,1) + t29 * mrSges(6,2) + t62 * t53 + t63 * t54 + t111;
t104 = -m(5) * t145 + t106;
t8 = t104 + (-t65 + t67) * t94 + (-mrSges(4,2) + mrSges(5,3)) * t93 + t127 * t76 + t132 * t49 + m(4) * t129;
t82 = (-mrSges(3,1) * t100 + mrSges(3,2) * t98) * qJD(1);
t4 = m(3) * (-t100 * g(3) - t134) - t83 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t82 * t125 + qJD(2) * t86 + t97 * t8 + t140 * t7;
t5 = m(3) * t119 - qJDD(2) * mrSges(3,2) + t84 * mrSges(3,3) - qJD(2) * t85 + t123 * t82 + t140 * t8 - t97 * t7;
t143 = t100 * t4 + t98 * t5;
t6 = m(2) * t120 + qJDD(1) * mrSges(2,1) - t102 * mrSges(2,2) - t146;
t1 = m(2) * t117 - t102 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t100 * t5 - t98 * t4;
t2 = [-m(1) * g(1) + t101 * t1 - t99 * t6, t1, t5, t8, -t49 * mrSges(5,2) - t76 * t66 - t113, t10, -t28 * mrSges(7,2) - t62 * t37 + t121; -m(1) * g(2) + t99 * t1 + t101 * t6, t6, t4, t7, t49 * mrSges(5,1) - t93 * mrSges(5,3) + t76 * t60 - t94 * t67 - t104, t11, t111; (-m(1) - m(2)) * g(3) + t143, -m(2) * g(3) + t143, t146, t105, t50 * mrSges(5,1) + t93 * mrSges(5,2) + t77 * t60 + t94 * t66 + t110, t106, -t47 * mrSges(7,1) + t29 * mrSges(7,2) + t63 * t37 - t74 * t52 + t141;];
f_new  = t2;
