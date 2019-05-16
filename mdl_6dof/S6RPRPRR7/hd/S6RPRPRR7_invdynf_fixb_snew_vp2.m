% Calculate vector of cutting forces with Newton-Euler
% S6RPRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2019-05-05 19:14
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRPRR7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR7_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR7_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR7_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR7_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR7_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 19:12:19
% EndTime: 2019-05-05 19:12:25
% DurationCPUTime: 2.51s
% Computational Cost: add. (33003->178), mult. (73315->233), div. (0->0), fcn. (51005->10), ass. (0->92)
t93 = sin(qJ(1));
t97 = cos(qJ(1));
t111 = -t97 * g(1) - t93 * g(2);
t108 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t111;
t127 = -m(2) - m(3);
t126 = -pkin(1) - pkin(7);
t125 = (mrSges(2,1) - mrSges(3,2));
t124 = -mrSges(2,2) + mrSges(3,3);
t119 = qJD(1) * qJD(3);
t92 = sin(qJ(3));
t113 = t92 * t119;
t115 = t93 * g(1) - t97 * g(2);
t98 = qJD(1) ^ 2;
t106 = -t98 * qJ(2) + qJDD(2) - t115;
t62 = t126 * qJDD(1) + t106;
t96 = cos(qJ(3));
t122 = t92 * g(3) + t96 * t62;
t76 = t96 * qJDD(1) - t113;
t38 = (-t76 - t113) * qJ(4) + (-t92 * t96 * t98 + qJDD(3)) * pkin(3) + t122;
t114 = -t96 * g(3) + t92 * t62;
t75 = -t92 * qJDD(1) - t96 * t119;
t120 = qJD(1) * t96;
t78 = qJD(3) * pkin(3) - qJ(4) * t120;
t87 = t92 ^ 2;
t39 = -t87 * t98 * pkin(3) + t75 * qJ(4) - qJD(3) * t78 + t114;
t88 = sin(pkin(10));
t89 = cos(pkin(10));
t69 = (-t88 * t92 + t89 * t96) * qJD(1);
t112 = -0.2e1 * qJD(4) * t69 + t89 * t38 - t88 * t39;
t54 = t88 * t75 + t89 * t76;
t68 = (-t88 * t96 - t89 * t92) * qJD(1);
t18 = (qJD(3) * t68 - t54) * pkin(8) + (t68 * t69 + qJDD(3)) * pkin(4) + t112;
t117 = 0.2e1 * qJD(4) * t68 + t88 * t38 + t89 * t39;
t53 = t89 * t75 - t88 * t76;
t61 = qJD(3) * pkin(4) - t69 * pkin(8);
t67 = t68 ^ 2;
t20 = -t67 * pkin(4) + t53 * pkin(8) - qJD(3) * t61 + t117;
t91 = sin(qJ(5));
t95 = cos(qJ(5));
t123 = t91 * t18 + t95 * t20;
t121 = qJD(1) * t92;
t51 = -t68 * mrSges(5,1) + t69 * mrSges(5,2);
t59 = -qJD(3) * mrSges(5,2) + t68 * mrSges(5,3);
t47 = t95 * t68 - t91 * t69;
t48 = t91 * t68 + t95 * t69;
t34 = -t47 * pkin(5) - t48 * pkin(9);
t85 = qJD(3) + qJD(5);
t83 = t85 ^ 2;
t84 = qJDD(3) + qJDD(5);
t15 = -t83 * pkin(5) + t84 * pkin(9) + t47 * t34 + t123;
t102 = -t75 * pkin(3) + qJDD(4) + t78 * t120 + (-qJ(4) * t87 + t126) * t98 + t108;
t101 = -t53 * pkin(4) - t67 * pkin(8) + t69 * t61 + t102;
t27 = -t48 * qJD(5) + t95 * t53 - t91 * t54;
t28 = t47 * qJD(5) + t91 * t53 + t95 * t54;
t16 = (-t47 * t85 - t28) * pkin(9) + (t48 * t85 - t27) * pkin(5) + t101;
t90 = sin(qJ(6));
t94 = cos(qJ(6));
t42 = -t90 * t48 + t94 * t85;
t22 = t42 * qJD(6) + t94 * t28 + t90 * t84;
t26 = qJDD(6) - t27;
t43 = t94 * t48 + t90 * t85;
t29 = -t42 * mrSges(7,1) + t43 * mrSges(7,2);
t46 = qJD(6) - t47;
t30 = -t46 * mrSges(7,2) + t42 * mrSges(7,3);
t12 = m(7) * (-t90 * t15 + t94 * t16) - t22 * mrSges(7,3) + t26 * mrSges(7,1) - t43 * t29 + t46 * t30;
t21 = -t43 * qJD(6) - t90 * t28 + t94 * t84;
t31 = t46 * mrSges(7,1) - t43 * mrSges(7,3);
t13 = m(7) * (t94 * t15 + t90 * t16) + t21 * mrSges(7,3) - t26 * mrSges(7,2) + t42 * t29 - t46 * t31;
t33 = -t47 * mrSges(6,1) + t48 * mrSges(6,2);
t45 = t85 * mrSges(6,1) - t48 * mrSges(6,3);
t8 = m(6) * t123 - t84 * mrSges(6,2) + t27 * mrSges(6,3) - t90 * t12 + t94 * t13 + t47 * t33 - t85 * t45;
t110 = t95 * t18 - t91 * t20;
t104 = m(7) * (-t84 * pkin(5) - t83 * pkin(9) + t48 * t34 - t110) - t21 * mrSges(7,1) + t22 * mrSges(7,2) - t42 * t30 + t43 * t31;
t44 = -t85 * mrSges(6,2) + t47 * mrSges(6,3);
t9 = m(6) * t110 + t84 * mrSges(6,1) - t28 * mrSges(6,3) - t48 * t33 + t85 * t44 - t104;
t5 = m(5) * t112 + qJDD(3) * mrSges(5,1) - t54 * mrSges(5,3) + qJD(3) * t59 - t69 * t51 + t91 * t8 + t95 * t9;
t60 = qJD(3) * mrSges(5,1) - t69 * mrSges(5,3);
t6 = m(5) * t117 - qJDD(3) * mrSges(5,2) + t53 * mrSges(5,3) - qJD(3) * t60 + t68 * t51 + t95 * t8 - t91 * t9;
t74 = (mrSges(4,1) * t92 + mrSges(4,2) * t96) * qJD(1);
t77 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t121;
t3 = m(4) * t122 + qJDD(3) * mrSges(4,1) - t76 * mrSges(4,3) + qJD(3) * t77 - t74 * t120 + t89 * t5 + t88 * t6;
t79 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t120;
t4 = m(4) * t114 - qJDD(3) * mrSges(4,2) + t75 * mrSges(4,3) - qJD(3) * t79 - t74 * t121 - t88 * t5 + t89 * t6;
t116 = -t92 * t3 + t96 * t4;
t107 = -m(3) * (-qJDD(1) * pkin(1) + t106) - t96 * t3 - t92 * t4;
t105 = m(6) * t101 - t27 * mrSges(6,1) + t28 * mrSges(6,2) + t94 * t12 + t90 * t13 - t47 * t44 + t48 * t45;
t103 = m(5) * t102 - t53 * mrSges(5,1) + t54 * mrSges(5,2) - t68 * t59 + t69 * t60 + t105;
t100 = -t75 * mrSges(4,1) + t103 + m(4) * (t126 * t98 + t108) + t77 * t121 + t79 * t120 + t76 * mrSges(4,2);
t99 = -m(3) * (t98 * pkin(1) - t108) + t100;
t7 = m(2) * t111 + t124 * qJDD(1) - (t125 * t98) + t99;
t1 = m(2) * t115 + t125 * qJDD(1) + t124 * t98 + t107;
t2 = [-m(1) * g(1) - t93 * t1 + t97 * t7, t7, -m(3) * g(3) + t116, t4, t6, t8, t13; -m(1) * g(2) + t97 * t1 + t93 * t7, t1, -(t98 * mrSges(3,2)) - qJDD(1) * mrSges(3,3) - t99, t3, t5, t9, t12; (-m(1) + t127) * g(3) + t116, t127 * g(3) + t116, qJDD(1) * mrSges(3,2) - t98 * mrSges(3,3) - t107, t100, t103, t105, t104;];
f_new  = t2;
