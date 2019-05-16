% Calculate vector of cutting forces with Newton-Euler
% S6RPRPRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-05-05 17:31
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRPRP1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP1_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP1_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP1_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP1_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP1_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:28:25
% EndTime: 2019-05-05 17:28:28
% DurationCPUTime: 1.83s
% Computational Cost: add. (21175->176), mult. (45514->225), div. (0->0), fcn. (29540->10), ass. (0->85)
t129 = -2 * qJD(4);
t92 = sin(qJ(1));
t95 = cos(qJ(1));
t111 = t92 * g(1) - t95 * g(2);
t72 = qJDD(1) * pkin(1) + t111;
t105 = -t95 * g(1) - t92 * g(2);
t97 = qJD(1) ^ 2;
t74 = -t97 * pkin(1) + t105;
t87 = sin(pkin(9));
t89 = cos(pkin(9));
t121 = t87 * t72 + t89 * t74;
t49 = -t97 * pkin(2) + qJDD(1) * pkin(7) + t121;
t85 = -g(3) + qJDD(2);
t91 = sin(qJ(3));
t94 = cos(qJ(3));
t107 = -t91 * t49 + t94 * t85;
t117 = qJD(1) * qJD(3);
t110 = t94 * t117;
t75 = t91 * qJDD(1) + t110;
t28 = (-t75 + t110) * qJ(4) + (t91 * t94 * t97 + qJDD(3)) * pkin(3) + t107;
t122 = t94 * t49 + t91 * t85;
t76 = t94 * qJDD(1) - t91 * t117;
t120 = qJD(1) * t91;
t77 = qJD(3) * pkin(3) - qJ(4) * t120;
t84 = t94 ^ 2;
t29 = -t84 * t97 * pkin(3) + t76 * qJ(4) - qJD(3) * t77 + t122;
t86 = sin(pkin(10));
t88 = cos(pkin(10));
t67 = (t86 * t94 + t88 * t91) * qJD(1);
t128 = t67 * t129 + t88 * t28 - t86 * t29;
t66 = (t86 * t91 - t88 * t94) * qJD(1);
t112 = t66 * t129 + t86 * t28 + t88 * t29;
t52 = t66 * pkin(4) - t67 * pkin(8);
t96 = qJD(3) ^ 2;
t21 = -t96 * pkin(4) + qJDD(3) * pkin(8) - t66 * t52 + t112;
t56 = -t86 * t75 + t88 * t76;
t57 = t88 * t75 + t86 * t76;
t106 = t89 * t72 - t87 * t74;
t102 = -qJDD(1) * pkin(2) - t106;
t99 = -t76 * pkin(3) + qJDD(4) + t77 * t120 + (-qJ(4) * t84 - pkin(7)) * t97 + t102;
t24 = (qJD(3) * t66 - t57) * pkin(8) + (qJD(3) * t67 - t56) * pkin(4) + t99;
t90 = sin(qJ(5));
t93 = cos(qJ(5));
t109 = -t90 * t21 + t93 * t24;
t59 = t93 * qJD(3) - t90 * t67;
t37 = t59 * qJD(5) + t90 * qJDD(3) + t93 * t57;
t65 = qJD(5) + t66;
t42 = -t65 * mrSges(7,2) + t59 * mrSges(7,3);
t55 = qJDD(5) - t56;
t60 = t90 * qJD(3) + t93 * t67;
t115 = m(7) * (-0.2e1 * qJD(6) * t60 + (t59 * t65 - t37) * qJ(6) + (t59 * t60 + t55) * pkin(5) + t109) + t65 * t42 + t55 * mrSges(7,1);
t39 = -t59 * mrSges(7,1) + t60 * mrSges(7,2);
t40 = -t59 * mrSges(6,1) + t60 * mrSges(6,2);
t43 = -t65 * mrSges(6,2) + t59 * mrSges(6,3);
t11 = m(6) * t109 + t55 * mrSges(6,1) + t65 * t43 + (-t40 - t39) * t60 + (-mrSges(6,3) - mrSges(7,3)) * t37 + t115;
t124 = t93 * t21 + t90 * t24;
t36 = -t60 * qJD(5) + t93 * qJDD(3) - t90 * t57;
t44 = t65 * pkin(5) - t60 * qJ(6);
t58 = t59 ^ 2;
t114 = m(7) * (-t58 * pkin(5) + t36 * qJ(6) + 0.2e1 * qJD(6) * t59 - t65 * t44 + t124) + t59 * t39 + t36 * mrSges(7,3);
t45 = t65 * mrSges(7,1) - t60 * mrSges(7,3);
t46 = t65 * mrSges(6,1) - t60 * mrSges(6,3);
t13 = m(6) * t124 + t36 * mrSges(6,3) + t59 * t40 + (-t46 - t45) * t65 + (-mrSges(6,2) - mrSges(7,2)) * t55 + t114;
t61 = -qJD(3) * mrSges(5,2) - t66 * mrSges(5,3);
t62 = qJD(3) * mrSges(5,1) - t67 * mrSges(5,3);
t101 = -m(5) * t99 + t56 * mrSges(5,1) - t57 * mrSges(5,2) - t93 * t11 - t90 * t13 - t66 * t61 - t67 * t62;
t78 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t120;
t119 = qJD(1) * t94;
t79 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t119;
t127 = (t91 * t78 - t94 * t79) * qJD(1) + m(4) * (-t97 * pkin(7) + t102) - t76 * mrSges(4,1) + t75 * mrSges(4,2) - t101;
t20 = -qJDD(3) * pkin(4) - t96 * pkin(8) + t67 * t52 - t128;
t113 = m(7) * (-t36 * pkin(5) - t58 * qJ(6) + t60 * t44 + qJDD(6) + t20) + t60 * t45 + t37 * mrSges(7,2);
t126 = -m(6) * t20 - t37 * mrSges(6,2) + (t43 + t42) * t59 + (mrSges(6,1) + mrSges(7,1)) * t36 - t60 * t46 - t113;
t51 = t66 * mrSges(5,1) + t67 * mrSges(5,2);
t14 = m(5) * t128 + qJDD(3) * mrSges(5,1) - t57 * mrSges(5,3) + qJD(3) * t61 - t67 * t51 + t126;
t73 = (-mrSges(4,1) * t94 + mrSges(4,2) * t91) * qJD(1);
t9 = m(5) * t112 - qJDD(3) * mrSges(5,2) + t56 * mrSges(5,3) - qJD(3) * t62 - t90 * t11 + t93 * t13 - t66 * t51;
t6 = m(4) * t107 + qJDD(3) * mrSges(4,1) - t75 * mrSges(4,3) + qJD(3) * t79 - t73 * t120 + t88 * t14 + t86 * t9;
t7 = m(4) * t122 - qJDD(3) * mrSges(4,2) + t76 * mrSges(4,3) - qJD(3) * t78 + t73 * t119 - t86 * t14 + t88 * t9;
t116 = m(3) * t85 + t94 * t6 + t91 * t7;
t8 = m(3) * t106 + qJDD(1) * mrSges(3,1) - t97 * mrSges(3,2) - t127;
t3 = m(3) * t121 - t97 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t91 * t6 + t94 * t7;
t2 = m(2) * t105 - t97 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t89 * t3 - t87 * t8;
t1 = m(2) * t111 + qJDD(1) * mrSges(2,1) - t97 * mrSges(2,2) + t87 * t3 + t89 * t8;
t4 = [-m(1) * g(1) - t92 * t1 + t95 * t2, t2, t3, t7, t9, t13, -t55 * mrSges(7,2) - t65 * t45 + t114; -m(1) * g(2) + t95 * t1 + t92 * t2, t1, t8, t6, t14, t11, -t37 * mrSges(7,3) - t60 * t39 + t115; (-m(1) - m(2)) * g(3) + t116, -m(2) * g(3) + t116, t116, t127, -t101, -t126, -t36 * mrSges(7,1) - t59 * t42 + t113;];
f_new  = t4;
