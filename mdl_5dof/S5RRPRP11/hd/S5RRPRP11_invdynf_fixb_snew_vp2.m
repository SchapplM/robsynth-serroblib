% Calculate vector of cutting forces with Newton-Euler
% S5RRPRP11
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
% Datum: 2019-12-31 20:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRPRP11_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP11_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP11_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP11_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP11_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP11_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP11_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP11_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP11_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:12:24
% EndTime: 2019-12-31 20:12:26
% DurationCPUTime: 0.69s
% Computational Cost: add. (4261->158), mult. (8632->191), div. (0->0), fcn. (4469->6), ass. (0->73)
t117 = -2 * qJD(3);
t73 = sin(qJ(2));
t99 = t73 * qJD(1);
t60 = mrSges(4,1) * t99 + qJD(2) * mrSges(4,2);
t101 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t99 - t60;
t108 = mrSges(3,1) - mrSges(4,2);
t79 = qJD(1) ^ 2;
t110 = t79 * pkin(6);
t76 = cos(qJ(2));
t98 = qJD(1) * qJD(2);
t92 = t76 * t98;
t54 = t73 * qJDD(1) + t92;
t93 = t73 * t98;
t55 = t76 * qJDD(1) - t93;
t100 = qJD(1) * t76;
t58 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t100;
t74 = sin(qJ(1));
t77 = cos(qJ(1));
t95 = t74 * g(1) - t77 * g(2);
t88 = -qJDD(1) * pkin(1) - t95;
t59 = -mrSges(4,1) * t100 - qJD(2) * mrSges(4,3);
t72 = sin(qJ(4));
t75 = cos(qJ(4));
t49 = t72 * qJD(2) + t75 * t100;
t50 = t75 * qJD(2) - t72 * t100;
t32 = t49 * mrSges(6,1) - t50 * mrSges(6,3);
t104 = -t49 * mrSges(5,1) - t50 * mrSges(5,2) - t32;
t61 = pkin(3) * t99 - qJD(2) * pkin(7);
t71 = t76 ^ 2;
t82 = pkin(2) * t93 + t99 * t117 + (-t54 - t92) * qJ(3) + t88;
t15 = -t61 * t99 + (-pkin(3) * t71 - pkin(6)) * t79 + (-pkin(2) - pkin(7)) * t55 + t82;
t91 = -t77 * g(1) - t74 * g(2);
t43 = -t79 * pkin(1) + qJDD(1) * pkin(6) + t91;
t103 = -t76 * g(3) - t73 * t43;
t51 = (-pkin(2) * t76 - qJ(3) * t73) * qJD(1);
t78 = qJD(2) ^ 2;
t23 = -qJDD(2) * pkin(2) - t78 * qJ(3) + t51 * t99 + qJDD(3) - t103;
t19 = (-t73 * t76 * t79 - qJDD(2)) * pkin(7) + (t54 - t92) * pkin(3) + t23;
t105 = t75 * t15 + t72 * t19;
t106 = -mrSges(5,3) - mrSges(6,2);
t28 = t50 * qJD(4) + t72 * qJDD(2) + t75 * t55;
t64 = qJD(4) + t99;
t36 = t64 * mrSges(5,1) - t50 * mrSges(5,3);
t48 = qJDD(4) + t54;
t31 = t49 * pkin(4) - t50 * qJ(5);
t37 = -t64 * mrSges(6,1) + t50 * mrSges(6,2);
t62 = t64 ^ 2;
t96 = m(6) * (-t62 * pkin(4) + t48 * qJ(5) + 0.2e1 * qJD(5) * t64 - t49 * t31 + t105) + t64 * t37 + t48 * mrSges(6,3);
t8 = m(5) * t105 - t48 * mrSges(5,2) + t104 * t49 + t106 * t28 - t64 * t36 + t96;
t90 = -t72 * t15 + t75 * t19;
t111 = m(6) * (-t48 * pkin(4) - t62 * qJ(5) + t50 * t31 + qJDD(5) - t90);
t29 = -t49 * qJD(4) + t75 * qJDD(2) - t72 * t55;
t34 = -t64 * mrSges(5,2) - t49 * mrSges(5,3);
t35 = -t49 * mrSges(6,2) + t64 * mrSges(6,3);
t9 = m(5) * t90 - t111 + (t34 + t35) * t64 + t104 * t50 + (mrSges(5,1) + mrSges(6,1)) * t48 + t106 * t29;
t89 = t72 * t9 - t75 * t8 - m(4) * (-t55 * pkin(2) - t110 + t82) - t59 * t100 + t54 * mrSges(4,3);
t116 = (t101 * t73 - t76 * t58) * qJD(1) - t108 * t55 + m(3) * (t88 - t110) + t54 * mrSges(3,2) - t89;
t94 = -t73 * g(3) + t76 * t43;
t115 = t78 * pkin(2) - qJDD(2) * qJ(3) + qJD(2) * t117 - t51 * t100 - t94;
t52 = (mrSges(4,2) * t76 - mrSges(4,3) * t73) * qJD(1);
t102 = t52 + (-mrSges(3,1) * t76 + mrSges(3,2) * t73) * qJD(1);
t107 = mrSges(3,3) + mrSges(4,1);
t87 = -m(4) * t23 - t72 * t8 - t75 * t9;
t4 = m(3) * t103 - t107 * t54 + t108 * qJDD(2) + (t58 - t59) * qJD(2) - t102 * t99 + t87;
t83 = -t71 * t79 * pkin(7) + t55 * pkin(3) + qJD(2) * t61 - t115;
t85 = -t29 * mrSges(6,3) - t50 * t37 + m(6) * (-0.2e1 * qJD(5) * t50 + (t49 * t64 - t29) * qJ(5) + (t50 * t64 + t28) * pkin(4) + t83) + t28 * mrSges(6,1) + t49 * t35;
t81 = m(5) * t83 + t28 * mrSges(5,1) + t29 * mrSges(5,2) + t49 * t34 + t50 * t36 + t85;
t80 = -m(4) * t115 + t81;
t6 = t80 + t107 * t55 - t101 * qJD(2) + (-mrSges(3,2) + mrSges(4,3)) * qJDD(2) + m(3) * t94 + t102 * t100;
t112 = t76 * t4 + t73 * t6;
t2 = m(2) * t95 + qJDD(1) * mrSges(2,1) - t79 * mrSges(2,2) - t116;
t1 = m(2) * t91 - t79 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t73 * t4 + t76 * t6;
t3 = [-m(1) * g(1) + t77 * t1 - t74 * t2, t1, t6, t55 * mrSges(4,2) - t60 * t99 - t89, t8, -t28 * mrSges(6,2) - t49 * t32 + t96; -m(1) * g(2) + t74 * t1 + t77 * t2, t2, t4, -t55 * mrSges(4,1) - qJDD(2) * mrSges(4,3) - qJD(2) * t60 - t52 * t100 - t80, t9, t85; (-m(1) - m(2)) * g(3) + t112, -m(2) * g(3) + t112, t116, t54 * mrSges(4,1) + qJDD(2) * mrSges(4,2) + qJD(2) * t59 + t52 * t99 - t87, t81, -t48 * mrSges(6,1) + t29 * mrSges(6,2) + t50 * t32 - t64 * t35 + t111;];
f_new = t3;
