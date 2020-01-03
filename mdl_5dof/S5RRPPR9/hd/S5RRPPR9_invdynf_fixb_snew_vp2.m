% Calculate vector of cutting forces with Newton-Euler
% S5RRPPR9
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% f_new [3x6]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRPPR9_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR9_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR9_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR9_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR9_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPPR9_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR9_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR9_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR9_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:40:43
% EndTime: 2019-12-31 19:40:44
% DurationCPUTime: 0.72s
% Computational Cost: add. (3244->167), mult. (6741->199), div. (0->0), fcn. (3126->6), ass. (0->78)
t112 = qJD(1) * qJD(2);
t90 = cos(qJ(2));
t107 = t90 * t112;
t88 = sin(qJ(1));
t91 = cos(qJ(1));
t116 = t88 * g(1) - t91 * g(2);
t93 = qJD(1) ^ 2;
t37 = -qJDD(1) * pkin(1) - t93 * pkin(6) - t116;
t87 = sin(qJ(2));
t57 = t87 * qJDD(1) + t107;
t108 = t87 * t112;
t58 = t90 * qJDD(1) - t108;
t103 = -t58 * pkin(2) + t37 + (-t107 - t57) * qJ(3);
t127 = pkin(3) + pkin(7);
t128 = -pkin(2) - pkin(7);
t113 = t87 * qJD(1);
t124 = t90 ^ 2 * t93;
t129 = 2 * qJD(3);
t62 = -qJD(2) * pkin(3) - qJ(4) * t113;
t95 = -qJ(4) * t124 + qJDD(4) - t103 + (t129 + t62) * t113;
t12 = t95 + t127 * t58 + (pkin(4) * t90 + t128 * t87) * t112 + t57 * pkin(4);
t105 = -t91 * g(1) - t88 * g(2);
t38 = -t93 * pkin(1) + qJDD(1) * pkin(6) + t105;
t120 = -t90 * g(3) - t87 * t38;
t52 = (-pkin(2) * t90 - qJ(3) * t87) * qJD(1);
t106 = t52 * t113 + qJDD(3) - t120;
t123 = t90 * t93;
t111 = qJD(1) * qJD(4);
t133 = -0.2e1 * t87 * t111 + (t107 - t57) * qJ(4);
t56 = (pkin(4) * t87 + pkin(7) * t90) * qJD(1);
t92 = qJD(2) ^ 2;
t15 = (-pkin(4) - qJ(3)) * t92 + (-pkin(3) * t123 - qJD(1) * t56) * t87 + (-pkin(2) - t127) * qJDD(2) + t106 + t133;
t114 = qJD(1) * t90;
t86 = sin(qJ(5));
t89 = cos(qJ(5));
t50 = -t89 * qJD(2) + t86 * t114;
t28 = t50 * qJD(5) - t86 * qJDD(2) - t89 * t58;
t51 = -t86 * qJD(2) - t89 * t114;
t29 = -t50 * mrSges(6,1) + t51 * mrSges(6,2);
t71 = qJD(5) + t113;
t30 = -t71 * mrSges(6,2) + t50 * mrSges(6,3);
t48 = qJDD(5) + t57;
t10 = m(6) * (t89 * t12 - t86 * t15) - t28 * mrSges(6,3) + t48 * mrSges(6,1) - t51 * t29 + t71 * t30;
t27 = -t51 * qJD(5) - t89 * qJDD(2) + t86 * t58;
t31 = t71 * mrSges(6,1) - t51 * mrSges(6,3);
t11 = m(6) * (t86 * t12 + t89 * t15) + t27 * mrSges(6,3) - t48 * mrSges(6,2) + t50 * t29 - t71 * t31;
t63 = qJD(2) * mrSges(5,2) - mrSges(5,3) * t113;
t66 = -qJD(2) * mrSges(5,1) + mrSges(5,3) * t114;
t104 = t89 * t10 + t86 * t11 + m(5) * (-pkin(2) * t108 + t58 * pkin(3) + t95) + t63 * t113 - t66 * t114 + t57 * mrSges(5,1) - t58 * mrSges(5,2);
t130 = -2 * qJD(3);
t100 = m(4) * ((pkin(2) * qJD(2) + t130) * t113 + t103) - t58 * mrSges(4,1) - t104;
t68 = mrSges(4,2) * t114 + qJD(2) * mrSges(4,3);
t117 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t114 + t68;
t64 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t113;
t65 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t113;
t135 = -(t117 * t90 - (t64 - t65) * t87) * qJD(1) + m(3) * t37 - t58 * mrSges(3,1) + (mrSges(3,2) - mrSges(4,3)) * t57 + t100;
t109 = -t87 * g(3) + t90 * t38;
t134 = qJDD(2) * qJ(3) + t52 * t114 + t109;
t55 = (mrSges(5,1) * t87 - mrSges(5,2) * t90) * qJD(1);
t119 = (-mrSges(3,1) * t90 + mrSges(3,2) * t87) * qJD(1) - t55;
t121 = mrSges(3,3) + mrSges(4,2);
t53 = (-mrSges(4,1) * t90 - mrSges(4,3) * t87) * qJD(1);
t24 = -qJDD(2) * pkin(2) - t92 * qJ(3) + t106;
t102 = t86 * t10 - t89 * t11 - m(5) * ((-t87 * t123 - qJDD(2)) * pkin(3) + t24 + t133) + t57 * mrSges(5,3) - qJD(2) * t66 - qJDD(2) * mrSges(5,2);
t99 = m(4) * t24 - t102;
t4 = m(3) * t120 - t121 * t57 + (mrSges(3,1) + mrSges(4,1)) * qJDD(2) + t117 * qJD(2) + (-t53 - t119) * t113 - t99;
t125 = t92 * pkin(2);
t75 = qJD(2) * t129;
t98 = pkin(3) * t124 + t58 * qJ(4) - t134;
t101 = -t27 * mrSges(6,1) - t50 * t30 + m(6) * (qJDD(2) * pkin(4) + qJD(2) * t62 + t75 + t128 * t92 + (-0.2e1 * qJD(4) - t56) * t114 - t98) + t28 * mrSges(6,2) + t51 * t31;
t97 = -m(5) * (0.2e1 * t90 * t111 + t125 + (t130 - t62) * qJD(2) + t98) + t101 - t58 * mrSges(5,3) + qJD(2) * t63 + qJDD(2) * mrSges(5,1);
t94 = m(4) * (t75 - t125 + t134) + t53 * t114 + qJD(2) * t65 + qJDD(2) * mrSges(4,3) + t97;
t6 = m(3) * t109 - qJDD(2) * mrSges(3,2) - qJD(2) * t64 + t119 * t114 + t121 * t58 + t94;
t126 = t90 * t4 + t87 * t6;
t110 = t55 * t114;
t2 = m(2) * t116 + qJDD(1) * mrSges(2,1) - t93 * mrSges(2,2) - t135;
t1 = m(2) * t105 - t93 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t87 * t4 + t90 * t6;
t3 = [-m(1) * g(1) + t91 * t1 - t88 * t2, t1, t6, t58 * mrSges(4,2) - t110 + t94, -t55 * t113 - t102, t11; -m(1) * g(2) + t88 * t1 + t91 * t2, t2, t4, -t57 * mrSges(4,3) + (-t87 * t65 - t90 * t68) * qJD(1) + t100, -t97 + t110, t10; (-m(1) - m(2)) * g(3) + t126, -m(2) * g(3) + t126, t135, -qJDD(2) * mrSges(4,1) + t57 * mrSges(4,2) - qJD(2) * t68 + (t53 - t55) * t113 + t99, t104, t101;];
f_new = t3;
