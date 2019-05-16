% Calculate vector of cutting forces with Newton-Euler
% S6RPPRRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
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
% Datum: 2019-05-05 15:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPPRRR5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR5_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR5_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR5_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR5_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR5_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:49:57
% EndTime: 2019-05-05 15:50:00
% DurationCPUTime: 1.02s
% Computational Cost: add. (10742->153), mult. (20807->190), div. (0->0), fcn. (12219->8), ass. (0->78)
t83 = sin(qJ(1));
t87 = cos(qJ(1));
t114 = t83 * g(1) - t87 * g(2);
t88 = qJD(1) ^ 2;
t51 = -qJDD(1) * pkin(1) - t88 * qJ(2) + qJDD(2) - t114;
t121 = -qJDD(1) * qJ(3) - (2 * qJD(3) * qJD(1)) + t51;
t102 = -t87 * g(1) - t83 * g(2);
t120 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t102;
t81 = sin(qJ(5));
t82 = sin(qJ(4));
t85 = cos(qJ(5));
t86 = cos(qJ(4));
t53 = (t81 * t86 + t82 * t85) * qJD(1);
t119 = -m(3) - m(4);
t118 = mrSges(3,2) - mrSges(4,3);
t117 = -mrSges(4,2) - mrSges(3,3);
t110 = qJD(1) * qJD(4);
t104 = t82 * t110;
t93 = qJDD(3) + (-pkin(1) - qJ(3)) * t88 + t120;
t41 = -qJDD(1) * pkin(7) + t93;
t115 = t82 * g(3) + t86 * t41;
t62 = t86 * qJDD(1) - t104;
t21 = (-t62 - t104) * pkin(8) + (-t82 * t86 * t88 + qJDD(4)) * pkin(4) + t115;
t105 = -t86 * g(3) + t82 * t41;
t61 = -t82 * qJDD(1) - t86 * t110;
t112 = qJD(1) * t86;
t65 = qJD(4) * pkin(4) - pkin(8) * t112;
t78 = t82 ^ 2;
t22 = -t78 * t88 * pkin(4) + t61 * pkin(8) - qJD(4) * t65 + t105;
t116 = t81 * t21 + t85 * t22;
t113 = qJD(1) * t82;
t111 = -m(2) + t119;
t107 = mrSges(2,1) - t118;
t60 = (mrSges(5,1) * t82 + mrSges(5,2) * t86) * qJD(1);
t63 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t113;
t54 = (-t81 * t82 + t85 * t86) * qJD(1);
t36 = t53 * pkin(5) - t54 * pkin(9);
t73 = qJD(4) + qJD(5);
t71 = t73 ^ 2;
t72 = qJDD(4) + qJDD(5);
t15 = -t71 * pkin(5) + t72 * pkin(9) - t53 * t36 + t116;
t29 = -t54 * qJD(5) + t85 * t61 - t81 * t62;
t30 = -t53 * qJD(5) + t81 * t61 + t85 * t62;
t89 = -t61 * pkin(4) + t65 * t112 + (-pkin(8) * t78 - pkin(7)) * t88 - t121;
t16 = (t53 * t73 - t30) * pkin(9) + (t54 * t73 - t29) * pkin(5) + t89;
t80 = sin(qJ(6));
t84 = cos(qJ(6));
t42 = -t80 * t54 + t84 * t73;
t18 = t42 * qJD(6) + t84 * t30 + t80 * t72;
t43 = t84 * t54 + t80 * t73;
t25 = -t42 * mrSges(7,1) + t43 * mrSges(7,2);
t28 = qJDD(6) - t29;
t52 = qJD(6) + t53;
t31 = -t52 * mrSges(7,2) + t42 * mrSges(7,3);
t12 = m(7) * (-t80 * t15 + t84 * t16) - t18 * mrSges(7,3) + t28 * mrSges(7,1) - t43 * t25 + t52 * t31;
t17 = -t43 * qJD(6) - t80 * t30 + t84 * t72;
t32 = t52 * mrSges(7,1) - t43 * mrSges(7,3);
t13 = m(7) * (t84 * t15 + t80 * t16) + t17 * mrSges(7,3) - t28 * mrSges(7,2) + t42 * t25 - t52 * t32;
t35 = t53 * mrSges(6,1) + t54 * mrSges(6,2);
t49 = t73 * mrSges(6,1) - t54 * mrSges(6,3);
t8 = m(6) * t116 - t72 * mrSges(6,2) + t29 * mrSges(6,3) - t80 * t12 + t84 * t13 - t53 * t35 - t73 * t49;
t101 = t85 * t21 - t81 * t22;
t48 = -t73 * mrSges(6,2) - t53 * mrSges(6,3);
t91 = m(7) * (-t72 * pkin(5) - t71 * pkin(9) + t54 * t36 - t101) - t17 * mrSges(7,1) + t18 * mrSges(7,2) - t42 * t31 + t43 * t32;
t9 = m(6) * t101 + t72 * mrSges(6,1) - t30 * mrSges(6,3) - t54 * t35 + t73 * t48 - t91;
t5 = m(5) * t115 + qJDD(4) * mrSges(5,1) - t62 * mrSges(5,3) + qJD(4) * t63 - t60 * t112 + t81 * t8 + t85 * t9;
t64 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t112;
t6 = m(5) * t105 - qJDD(4) * mrSges(5,2) + t61 * mrSges(5,3) - qJD(4) * t64 - t60 * t113 + t85 * t8 - t81 * t9;
t106 = -t82 * t5 + t86 * t6;
t103 = m(4) * t93 + qJDD(1) * mrSges(4,2) + t86 * t5 + t82 * t6;
t98 = m(3) * (t88 * pkin(1) - t120) - t103;
t96 = m(6) * t89 - t29 * mrSges(6,1) + t30 * mrSges(6,2) + t84 * t12 + t80 * t13 + t53 * t48 + t54 * t49;
t94 = m(5) * (-t88 * pkin(7) - t121) + t62 * mrSges(5,2) - t61 * mrSges(5,1) + t64 * t112 + t63 * t113 + t96;
t92 = m(4) * t121 - t94;
t90 = m(3) * t51 + t92;
t7 = (-mrSges(2,2) - t117) * t88 + t107 * qJDD(1) + m(2) * t114 - t90;
t1 = m(2) * t102 + (-mrSges(2,2) + mrSges(3,3)) * qJDD(1) - t107 * t88 - t98;
t2 = [-m(1) * g(1) + t87 * t1 - t83 * t7, t1, t119 * g(3) + t106, -m(4) * g(3) + t106, t6, t8, t13; -m(1) * g(2) + t83 * t1 + t87 * t7, t7, -qJDD(1) * mrSges(3,3) - t118 * t88 + t98, -t88 * mrSges(4,2) - qJDD(1) * mrSges(4,3) + t92, t5, t9, t12; (-m(1) + t111) * g(3) + t106, t111 * g(3) + t106, t118 * qJDD(1) + t117 * t88 + t90, -t88 * mrSges(4,3) + t103, t94, t96, t91;];
f_new  = t2;
