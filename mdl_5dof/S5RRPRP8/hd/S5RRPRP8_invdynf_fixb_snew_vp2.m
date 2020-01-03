% Calculate vector of cutting forces with Newton-Euler
% S5RRPRP8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRPRP8_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP8_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP8_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP8_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP8_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP8_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP8_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:03:22
% EndTime: 2019-12-31 20:03:24
% DurationCPUTime: 0.72s
% Computational Cost: add. (4600->157), mult. (9653->194), div. (0->0), fcn. (5215->6), ass. (0->70)
t89 = cos(qJ(2));
t110 = qJD(1) * t89;
t68 = mrSges(4,2) * t110 + qJD(2) * mrSges(4,3);
t113 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t110 + t68;
t87 = sin(qJ(1));
t90 = cos(qJ(1));
t112 = t87 * g(1) - t90 * g(2);
t92 = qJD(1) ^ 2;
t53 = -qJDD(1) * pkin(1) - t92 * pkin(6) - t112;
t109 = qJD(1) * qJD(2);
t104 = t89 * t109;
t86 = sin(qJ(2));
t62 = t86 * qJDD(1) + t104;
t105 = t86 * t109;
t63 = t89 * qJDD(1) - t105;
t111 = qJD(1) * t86;
t65 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t111;
t66 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t111;
t100 = -t63 * pkin(2) + t53 + (-t104 - t62) * qJ(3);
t85 = sin(qJ(4));
t88 = cos(qJ(4));
t52 = (-t85 * t89 + t86 * t88) * qJD(1);
t31 = -t52 * qJD(4) - t85 * t62 - t88 * t63;
t51 = (-t85 * t86 - t88 * t89) * qJD(1);
t32 = t51 * qJD(4) + t88 * t62 - t85 * t63;
t79 = -qJD(2) + qJD(4);
t40 = -t79 * mrSges(6,2) + t51 * mrSges(6,3);
t42 = t79 * pkin(4) - t52 * qJ(5);
t43 = t79 * mrSges(6,1) - t52 * mrSges(6,3);
t50 = t51 ^ 2;
t119 = t89 ^ 2 * t92;
t121 = 2 * qJD(3);
t69 = -qJD(2) * pkin(3) - pkin(7) * t111;
t93 = -pkin(2) * t105 + t63 * pkin(3) - pkin(7) * t119 - t100 + (t121 + t69) * t111;
t102 = m(6) * (-t31 * pkin(4) - t50 * qJ(5) + t52 * t42 + qJDD(5) + t93) + t32 * mrSges(6,2) - t31 * mrSges(6,1) + t52 * t43 - t51 * t40;
t41 = -t79 * mrSges(5,2) + t51 * mrSges(5,3);
t44 = t79 * mrSges(5,1) - t52 * mrSges(5,3);
t96 = m(5) * t93 - t31 * mrSges(5,1) + t32 * mrSges(5,2) - t51 * t41 + t52 * t44 + t102;
t95 = m(4) * ((pkin(2) * qJD(2) - (2 * qJD(3))) * t111 + t100) - t63 * mrSges(4,1) - t96;
t124 = -(t113 * t89 - (t65 - t66) * t86) * qJD(1) + m(3) * t53 - t63 * mrSges(3,1) + (mrSges(3,2) - mrSges(4,3)) * t62 + t95;
t101 = -t90 * g(1) - t87 * g(2);
t54 = -t92 * pkin(1) + qJDD(1) * pkin(6) + t101;
t106 = -t86 * g(3) + t89 * t54;
t117 = mrSges(3,3) + mrSges(4,2);
t61 = (-mrSges(3,1) * t89 + mrSges(3,2) * t86) * qJD(1);
t60 = (-mrSges(4,1) * t89 - mrSges(4,3) * t86) * qJD(1);
t59 = (-pkin(2) * t89 - qJ(3) * t86) * qJD(1);
t91 = qJD(2) ^ 2;
t97 = -t91 * pkin(2) + qJDD(2) * qJ(3) + qJD(2) * t121 + t59 * t110 + t106;
t19 = -pkin(3) * t119 - t63 * pkin(7) + qJD(2) * t69 + t97;
t115 = -t89 * g(3) - t86 * t54;
t30 = -qJDD(2) * pkin(2) - t91 * qJ(3) + t59 * t111 + qJDD(3) - t115;
t20 = (-t62 + t104) * pkin(7) + (-t86 * t89 * t92 - qJDD(2)) * pkin(3) + t30;
t103 = -t85 * t19 + t88 * t20;
t78 = -qJDD(2) + qJDD(4);
t108 = m(6) * (-0.2e1 * qJD(5) * t52 + (t51 * t79 - t32) * qJ(5) + (t51 * t52 + t78) * pkin(4) + t103) + t79 * t40 + t78 * mrSges(6,1);
t37 = -t51 * mrSges(6,1) + t52 * mrSges(6,2);
t38 = -t51 * mrSges(5,1) + t52 * mrSges(5,2);
t7 = m(5) * t103 + t78 * mrSges(5,1) + t79 * t41 + (-t38 - t37) * t52 + (-mrSges(5,3) - mrSges(6,3)) * t32 + t108;
t116 = t88 * t19 + t85 * t20;
t107 = m(6) * (-t50 * pkin(4) + t31 * qJ(5) + 0.2e1 * qJD(5) * t51 - t79 * t42 + t116) + t51 * t37 + t31 * mrSges(6,3);
t9 = m(5) * t116 + t31 * mrSges(5,3) + t51 * t38 + (-t44 - t43) * t79 + (-mrSges(5,2) - mrSges(6,2)) * t78 + t107;
t99 = m(4) * t97 + qJDD(2) * mrSges(4,3) + qJD(2) * t66 + t60 * t110 - t85 * t7 + t88 * t9;
t4 = m(3) * t106 - qJDD(2) * mrSges(3,2) - qJD(2) * t65 + t61 * t110 + t117 * t63 + t99;
t98 = -m(4) * t30 - t88 * t7 - t85 * t9;
t5 = m(3) * t115 - t117 * t62 + (mrSges(3,1) + mrSges(4,1)) * qJDD(2) + t113 * qJD(2) + (-t60 - t61) * t111 + t98;
t120 = t86 * t4 + t89 * t5;
t6 = m(2) * t112 + qJDD(1) * mrSges(2,1) - t92 * mrSges(2,2) - t124;
t1 = m(2) * t101 - t92 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t89 * t4 - t86 * t5;
t2 = [-m(1) * g(1) + t90 * t1 - t87 * t6, t1, t4, t63 * mrSges(4,2) + t99, t9, -t78 * mrSges(6,2) - t79 * t43 + t107; -m(1) * g(2) + t87 * t1 + t90 * t6, t6, t5, t95 + (-t86 * t66 - t89 * t68) * qJD(1) - t62 * mrSges(4,3), t7, -t32 * mrSges(6,3) - t52 * t37 + t108; (-m(1) - m(2)) * g(3) + t120, -m(2) * g(3) + t120, t124, -qJDD(2) * mrSges(4,1) + t62 * mrSges(4,2) - qJD(2) * t68 + t60 * t111 - t98, t96, t102;];
f_new = t2;
