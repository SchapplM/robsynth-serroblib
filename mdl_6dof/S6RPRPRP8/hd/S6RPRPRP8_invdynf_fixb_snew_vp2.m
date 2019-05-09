% Calculate vector of cutting forces with Newton-Euler
% S6RPRPRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2019-05-05 18:03
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRPRP8_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP8_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP8_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP8_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP8_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP8_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP8_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:00:28
% EndTime: 2019-05-05 18:00:31
% DurationCPUTime: 1.41s
% Computational Cost: add. (14139->175), mult. (30563->217), div. (0->0), fcn. (19483->8), ass. (0->84)
t132 = -2 * qJD(4);
t88 = sin(qJ(1));
t90 = cos(qJ(1));
t104 = -t90 * g(1) - t88 * g(2);
t100 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t104;
t114 = qJD(1) * qJD(3);
t87 = sin(qJ(3));
t106 = t87 * t114;
t128 = -pkin(1) - pkin(7);
t108 = t88 * g(1) - t90 * g(2);
t92 = qJD(1) ^ 2;
t98 = -t92 * qJ(2) + qJDD(2) - t108;
t57 = t128 * qJDD(1) + t98;
t89 = cos(qJ(3));
t118 = t87 * g(3) + t89 * t57;
t72 = t89 * qJDD(1) - t106;
t29 = (-t72 - t106) * qJ(4) + (-t87 * t89 * t92 + qJDD(3)) * pkin(3) + t118;
t107 = -t89 * g(3) + t87 * t57;
t71 = -t87 * qJDD(1) - t114 * t89;
t116 = qJD(1) * t89;
t74 = qJD(3) * pkin(3) - qJ(4) * t116;
t83 = t87 ^ 2;
t30 = -t83 * t92 * pkin(3) + t71 * qJ(4) - qJD(3) * t74 + t107;
t117 = qJD(1) * t87;
t84 = sin(pkin(9));
t85 = cos(pkin(9));
t65 = t116 * t85 - t117 * t84;
t131 = t65 * t132 + t85 * t29 - t84 * t30;
t64 = (t84 * t89 + t85 * t87) * qJD(1);
t43 = t64 * pkin(4) - t65 * pkin(8);
t91 = qJD(3) ^ 2;
t17 = -qJDD(3) * pkin(4) - t91 * pkin(8) + t65 * t43 - t131;
t126 = cos(qJ(5));
t48 = t84 * t71 + t85 * t72;
t86 = sin(qJ(5));
t50 = t86 * qJD(3) + t126 * t65;
t23 = t50 * qJD(5) - t126 * qJDD(3) + t86 * t48;
t49 = -t126 * qJD(3) + t86 * t65;
t24 = -t49 * qJD(5) + t86 * qJDD(3) + t126 * t48;
t62 = qJD(5) + t64;
t37 = -t49 * mrSges(7,2) + t62 * mrSges(7,3);
t111 = m(7) * (-0.2e1 * qJD(6) * t50 + (t49 * t62 - t24) * qJ(6) + (t50 * t62 + t23) * pkin(5) + t17) + t49 * t37 + t23 * mrSges(7,1);
t38 = -t62 * mrSges(6,2) - t49 * mrSges(6,3);
t39 = t62 * mrSges(6,1) - t50 * mrSges(6,3);
t40 = -t62 * mrSges(7,1) + t50 * mrSges(7,2);
t130 = m(6) * t17 + t23 * mrSges(6,1) + (t39 - t40) * t50 + (mrSges(6,2) - mrSges(7,3)) * t24 + t49 * t38 + t111;
t129 = -m(2) - m(3);
t110 = t64 * t132 + t84 * t29 + t85 * t30;
t18 = -t91 * pkin(4) + qJDD(3) * pkin(8) - t64 * t43 + t110;
t47 = t85 * t71 - t84 * t72;
t94 = -t71 * pkin(3) + qJDD(4) + t74 * t116 + (-qJ(4) * t83 + t128) * t92 + t100;
t20 = (qJD(3) * t64 - t48) * pkin(8) + (qJD(3) * t65 - t47) * pkin(4) + t94;
t101 = t126 * t20 - t86 * t18;
t33 = t49 * pkin(5) - t50 * qJ(6);
t46 = qJDD(5) - t47;
t61 = t62 ^ 2;
t127 = m(7) * (-t46 * pkin(5) - t61 * qJ(6) + t50 * t33 + qJDD(6) - t101);
t125 = (mrSges(2,1) - mrSges(3,2));
t124 = -mrSges(2,2) + mrSges(3,3);
t122 = -mrSges(6,3) - mrSges(7,2);
t121 = t126 * t18 + t86 * t20;
t34 = t49 * mrSges(7,1) - t50 * mrSges(7,3);
t120 = -t49 * mrSges(6,1) - t50 * mrSges(6,2) - t34;
t112 = m(7) * (-t61 * pkin(5) + t46 * qJ(6) + 0.2e1 * qJD(6) * t62 - t49 * t33 + t121) + t62 * t40 + t46 * mrSges(7,3);
t11 = m(6) * t101 - t127 + (t38 + t37) * t62 + t120 * t50 + (mrSges(6,1) + mrSges(7,1)) * t46 + t122 * t24;
t42 = t64 * mrSges(5,1) + t65 * mrSges(5,2);
t56 = qJD(3) * mrSges(5,1) - t65 * mrSges(5,3);
t9 = m(6) * t121 - t46 * mrSges(6,2) + t120 * t49 + t122 * t23 - t62 * t39 + t112;
t6 = m(5) * t110 - qJDD(3) * mrSges(5,2) + t47 * mrSges(5,3) - qJD(3) * t56 - t86 * t11 + t126 * t9 - t64 * t42;
t55 = -qJD(3) * mrSges(5,2) - t64 * mrSges(5,3);
t7 = m(5) * t131 + qJDD(3) * mrSges(5,1) - t48 * mrSges(5,3) + qJD(3) * t55 - t65 * t42 - t130;
t70 = (mrSges(4,1) * t87 + mrSges(4,2) * t89) * qJD(1);
t73 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t117;
t3 = m(4) * t118 + qJDD(3) * mrSges(4,1) - t72 * mrSges(4,3) + qJD(3) * t73 - t116 * t70 + t84 * t6 + t85 * t7;
t75 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t116;
t4 = m(4) * t107 - qJDD(3) * mrSges(4,2) + t71 * mrSges(4,3) - qJD(3) * t75 - t117 * t70 + t85 * t6 - t84 * t7;
t109 = -t87 * t3 + t89 * t4;
t99 = -m(3) * (-qJDD(1) * pkin(1) + t98) - t89 * t3 - t87 * t4;
t97 = m(5) * t94 - t47 * mrSges(5,1) + t48 * mrSges(5,2) + t126 * t11 + t64 * t55 + t65 * t56 + t86 * t9;
t95 = -t71 * mrSges(4,1) + m(4) * (t128 * t92 + t100) + t73 * t117 + t75 * t116 + t72 * mrSges(4,2) + t97;
t93 = -m(3) * (t92 * pkin(1) - t100) + t95;
t5 = m(2) * t104 + t124 * qJDD(1) - (t125 * t92) + t93;
t1 = m(2) * t108 + t125 * qJDD(1) + t124 * t92 + t99;
t2 = [-m(1) * g(1) - t88 * t1 + t90 * t5, t5, -m(3) * g(3) + t109, t4, t6, t9, -t23 * mrSges(7,2) - t49 * t34 + t112; -m(1) * g(2) + t90 * t1 + t88 * t5, t1, -(t92 * mrSges(3,2)) - qJDD(1) * mrSges(3,3) - t93, t3, t7, t11, -t24 * mrSges(7,3) - t50 * t40 + t111; (-m(1) + t129) * g(3) + t109, t129 * g(3) + t109, qJDD(1) * mrSges(3,2) - t92 * mrSges(3,3) - t99, t95, t97, t130, -t46 * mrSges(7,1) + t24 * mrSges(7,2) + t50 * t34 - t62 * t37 + t127;];
f_new  = t2;
