% Calculate vector of cutting forces with Newton-Euler
% S6RPRPRP3
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
% Datum: 2019-05-05 17:41
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRPRP3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP3_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP3_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP3_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP3_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP3_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:37:27
% EndTime: 2019-05-05 17:37:31
% DurationCPUTime: 1.88s
% Computational Cost: add. (24410->172), mult. (50153->221), div. (0->0), fcn. (31969->10), ass. (0->86)
t90 = sin(qJ(1));
t92 = cos(qJ(1));
t107 = t90 * g(1) - t92 * g(2);
t69 = qJDD(1) * pkin(1) + t107;
t102 = -t92 * g(1) - t90 * g(2);
t94 = qJD(1) ^ 2;
t72 = -t94 * pkin(1) + t102;
t85 = sin(pkin(9));
t87 = cos(pkin(9));
t104 = t87 * t69 - t85 * t72;
t45 = -qJDD(1) * pkin(2) - t94 * pkin(7) - t104;
t119 = cos(qJ(5));
t111 = qJD(1) * qJD(3);
t91 = cos(qJ(3));
t106 = t91 * t111;
t89 = sin(qJ(3));
t73 = t89 * qJDD(1) + t106;
t80 = t89 * t111;
t74 = t91 * qJDD(1) - t80;
t31 = (-t73 - t106) * qJ(4) + (-t74 + t80) * pkin(3) + t45;
t112 = t91 * qJD(1);
t114 = t85 * t69 + t87 * t72;
t46 = -t94 * pkin(2) + qJDD(1) * pkin(7) + t114;
t83 = -g(3) + qJDD(2);
t115 = t91 * t46 + t89 * t83;
t70 = (-pkin(3) * t91 - qJ(4) * t89) * qJD(1);
t93 = qJD(3) ^ 2;
t37 = -t93 * pkin(3) + qJDD(3) * qJ(4) + t70 * t112 + t115;
t113 = qJD(1) * t89;
t84 = sin(pkin(10));
t86 = cos(pkin(10));
t67 = t84 * qJD(3) + t86 * t113;
t103 = -0.2e1 * qJD(4) * t67 + t86 * t31 - t84 * t37;
t55 = t84 * qJDD(3) + t86 * t73;
t66 = t86 * qJD(3) - t84 * t113;
t18 = (-t66 * t112 - t55) * pkin(8) + (t66 * t67 - t74) * pkin(4) + t103;
t108 = 0.2e1 * qJD(4) * t66 + t84 * t31 + t86 * t37;
t54 = t86 * qJDD(3) - t84 * t73;
t56 = -pkin(4) * t112 - t67 * pkin(8);
t65 = t66 ^ 2;
t20 = -t65 * pkin(4) + t54 * pkin(8) + t56 * t112 + t108;
t88 = sin(qJ(5));
t100 = t119 * t18 - t88 * t20;
t48 = -t119 * t66 + t88 * t67;
t49 = t119 * t67 + t88 * t66;
t35 = t48 * mrSges(7,1) - t49 * mrSges(7,3);
t116 = -t48 * mrSges(6,1) - t49 * mrSges(6,2) - t35;
t118 = -mrSges(6,3) - mrSges(7,2);
t34 = t48 * pkin(5) - t49 * qJ(6);
t68 = qJDD(5) - t74;
t78 = qJD(5) - t112;
t77 = t78 ^ 2;
t120 = m(7) * (-t68 * pkin(5) - t77 * qJ(6) + t49 * t34 + qJDD(6) - t100);
t26 = -t48 * qJD(5) + t119 * t55 + t88 * t54;
t39 = -t78 * mrSges(6,2) - t48 * mrSges(6,3);
t42 = -t48 * mrSges(7,2) + t78 * mrSges(7,3);
t10 = m(6) * t100 - t120 + (t39 + t42) * t78 + (mrSges(6,1) + mrSges(7,1)) * t68 + t116 * t49 + t118 * t26;
t50 = -t66 * mrSges(5,1) + t67 * mrSges(5,2);
t52 = mrSges(5,2) * t112 + t66 * mrSges(5,3);
t117 = t119 * t20 + t88 * t18;
t41 = -t78 * mrSges(7,1) + t49 * mrSges(7,2);
t109 = m(7) * (-t77 * pkin(5) + t68 * qJ(6) + 0.2e1 * qJD(6) * t78 - t48 * t34 + t117) + t78 * t41 + t68 * mrSges(7,3);
t25 = t49 * qJD(5) - t119 * t54 + t88 * t55;
t40 = t78 * mrSges(6,1) - t49 * mrSges(6,3);
t9 = m(6) * t117 - t68 * mrSges(6,2) + t116 * t48 + t118 * t25 - t78 * t40 + t109;
t7 = m(5) * t103 - t74 * mrSges(5,1) - t55 * mrSges(5,3) + t119 * t10 - t52 * t112 - t67 * t50 + t88 * t9;
t75 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t113;
t76 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t112;
t53 = -mrSges(5,1) * t112 - t67 * mrSges(5,3);
t8 = m(5) * t108 + t74 * mrSges(5,2) + t54 * mrSges(5,3) - t88 * t10 + t53 * t112 + t119 * t9 + t66 * t50;
t121 = m(4) * t45 - t74 * mrSges(4,1) + t73 * mrSges(4,2) + t86 * t7 + t84 * t8 + (t89 * t75 - t91 * t76) * qJD(1);
t105 = -t89 * t46 + t91 * t83;
t71 = (-mrSges(4,1) * t91 + mrSges(4,2) * t89) * qJD(1);
t33 = -qJDD(3) * pkin(3) - t93 * qJ(4) + t70 * t113 + qJDD(4) - t105;
t96 = -t54 * pkin(4) - t65 * pkin(8) + t67 * t56 + t33;
t99 = t26 * mrSges(7,3) + t49 * t41 - m(7) * (-0.2e1 * qJD(6) * t49 + (t48 * t78 - t26) * qJ(6) + (t49 * t78 + t25) * pkin(5) + t96) - t25 * mrSges(7,1) - t48 * t42;
t97 = m(6) * t96 + t25 * mrSges(6,1) + t26 * mrSges(6,2) + t48 * t39 + t49 * t40 - t99;
t95 = m(5) * t33 - t54 * mrSges(5,1) + t55 * mrSges(5,2) - t66 * t52 + t67 * t53 + t97;
t12 = m(4) * t105 + qJDD(3) * mrSges(4,1) - t73 * mrSges(4,3) + qJD(3) * t76 - t71 * t113 - t95;
t6 = m(4) * t115 - qJDD(3) * mrSges(4,2) + t74 * mrSges(4,3) - qJD(3) * t75 + t71 * t112 - t84 * t7 + t86 * t8;
t110 = m(3) * t83 + t91 * t12 + t89 * t6;
t4 = m(3) * t104 + qJDD(1) * mrSges(3,1) - t94 * mrSges(3,2) - t121;
t3 = m(3) * t114 - t94 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t89 * t12 + t91 * t6;
t2 = m(2) * t102 - t94 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t87 * t3 - t85 * t4;
t1 = m(2) * t107 + qJDD(1) * mrSges(2,1) - t94 * mrSges(2,2) + t85 * t3 + t87 * t4;
t5 = [-m(1) * g(1) - t90 * t1 + t92 * t2, t2, t3, t6, t8, t9, -t25 * mrSges(7,2) - t48 * t35 + t109; -m(1) * g(2) + t92 * t1 + t90 * t2, t1, t4, t12, t7, t10, -t99; (-m(1) - m(2)) * g(3) + t110, -m(2) * g(3) + t110, t110, t121, t95, t97, -t68 * mrSges(7,1) + t26 * mrSges(7,2) + t49 * t35 - t78 * t42 + t120;];
f_new  = t5;
