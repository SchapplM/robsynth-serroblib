% Calculate vector of cutting forces with Newton-Euler
% S6PRRRPP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-05-05 06:48
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRRRPP2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP2_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPP2_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP2_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP2_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPP2_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 06:42:03
% EndTime: 2019-05-05 06:42:07
% DurationCPUTime: 1.25s
% Computational Cost: add. (14036->174), mult. (26544->210), div. (0->0), fcn. (17218->10), ass. (0->86)
t130 = cos(qJ(4));
t96 = cos(qJ(3));
t118 = t96 * qJD(2);
t90 = sin(pkin(6));
t95 = sin(qJ(2));
t127 = t90 * t95;
t89 = sin(pkin(10));
t91 = cos(pkin(10));
t75 = t89 * g(1) - t91 * g(2);
t92 = cos(pkin(6));
t128 = t75 * t92;
t76 = -t91 * g(1) - t89 * g(2);
t88 = -g(3) + qJDD(1);
t97 = cos(qJ(2));
t114 = t88 * t127 + t95 * t128 + t97 * t76;
t99 = qJD(2) ^ 2;
t31 = -t99 * pkin(2) + qJDD(2) * pkin(8) + t114;
t59 = -t90 * t75 + t92 * t88;
t94 = sin(qJ(3));
t122 = t96 * t31 + t94 * t59;
t72 = (-pkin(3) * t96 - pkin(9) * t94) * qJD(2);
t98 = qJD(3) ^ 2;
t25 = -t98 * pkin(3) + qJDD(3) * pkin(9) + t72 * t118 + t122;
t117 = qJD(2) * qJD(3);
t112 = t96 * t117;
t136 = (t88 * t90 + t128) * t97 - t95 * t76;
t30 = -qJDD(2) * pkin(2) - t99 * pkin(8) - t136;
t73 = t94 * qJDD(2) + t112;
t84 = t94 * t117;
t74 = t96 * qJDD(2) - t84;
t27 = (-t73 - t112) * pkin(9) + (-t74 + t84) * pkin(3) + t30;
t93 = sin(qJ(4));
t124 = t130 * t25 + t93 * t27;
t132 = -2 * qJD(5);
t119 = qJD(2) * t94;
t69 = -t130 * qJD(3) + t93 * t119;
t70 = t93 * qJD(3) + t130 * t119;
t46 = t69 * pkin(4) - t70 * qJ(5);
t66 = -qJDD(4) + t74;
t82 = -qJD(4) + t118;
t81 = t82 ^ 2;
t105 = -t81 * pkin(4) - t66 * qJ(5) + t82 * t132 - t69 * t46 + t124;
t57 = t82 * mrSges(6,1) + t70 * mrSges(6,2);
t137 = m(6) * t105 - t66 * mrSges(6,3) - t82 * t57;
t123 = -t94 * t31 + t96 * t59;
t102 = qJDD(3) * pkin(3) + t98 * pkin(9) - t72 * t119 + t123;
t129 = t69 * t82;
t43 = -t69 * qJD(4) + t93 * qJDD(3) + t130 * t73;
t135 = -(t43 + t129) * qJ(5) - t102;
t109 = t130 * t27 - t93 * t25;
t19 = t66 * pkin(4) - t81 * qJ(5) + t70 * t46 + qJDD(5) - t109;
t52 = -t82 * mrSges(7,2) + t69 * mrSges(7,3);
t116 = m(7) * (-0.2e1 * qJD(6) * t70 + (-t43 + t129) * qJ(6) + (t69 * t70 + t66) * pkin(5) + t19) + t82 * t52 + t66 * mrSges(7,1);
t106 = m(6) * t19 + t116;
t47 = t69 * mrSges(6,1) - t70 * mrSges(6,3);
t121 = -t69 * mrSges(5,1) - t70 * mrSges(5,2) - t47;
t125 = -mrSges(5,3) - mrSges(6,2);
t48 = -t69 * mrSges(7,1) + t70 * mrSges(7,2);
t53 = t82 * mrSges(5,2) - t69 * mrSges(5,3);
t58 = -t69 * mrSges(6,2) - t82 * mrSges(6,3);
t11 = m(5) * t109 + (-t53 - t58) * t82 + (-mrSges(5,1) - mrSges(6,1)) * t66 + (t48 + t121) * t70 + (mrSges(7,3) + t125) * t43 - t106;
t42 = t70 * qJD(4) - t130 * qJDD(3) + t93 * t73;
t54 = t82 * pkin(5) - t70 * qJ(6);
t65 = t69 ^ 2;
t115 = m(7) * (-t65 * pkin(5) + t42 * qJ(6) + 0.2e1 * qJD(6) * t69 - t82 * t54 + t105) + t69 * t48 + t42 * mrSges(7,3);
t55 = t82 * mrSges(7,1) - t70 * mrSges(7,3);
t56 = -t82 * mrSges(5,1) - t70 * mrSges(5,3);
t12 = m(5) * t124 + (t56 - t55) * t82 + t121 * t69 + (mrSges(5,2) - mrSges(7,2)) * t66 + t125 * t42 + t115 + t137;
t77 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t119;
t78 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t118;
t134 = m(4) * t30 - t74 * mrSges(4,1) + t73 * mrSges(4,2) + (t77 * t94 - t78 * t96) * qJD(2) + t130 * t11 + t93 * t12;
t110 = m(7) * (-t65 * qJ(6) + qJDD(6) + (-pkin(4) - pkin(5)) * t42 + (pkin(4) * t82 + (2 * qJD(5)) + t54) * t70 - t135) + t43 * mrSges(7,2) - t42 * mrSges(7,1) + t70 * t55 - t69 * t52;
t104 = m(6) * (t70 * t132 + (-t70 * t82 + t42) * pkin(4) + t135) + t42 * mrSges(6,1) + t69 * t58 - t110;
t133 = -m(5) * t102 + t42 * mrSges(5,1) + (t56 - t57) * t70 + (mrSges(5,2) - mrSges(6,3)) * t43 + t69 * t53 + t104;
t8 = m(3) * t136 + qJDD(2) * mrSges(3,1) - t99 * mrSges(3,2) - t134;
t131 = t8 * t97;
t71 = (-mrSges(4,1) * t96 + mrSges(4,2) * t94) * qJD(2);
t10 = m(4) * t123 + qJDD(3) * mrSges(4,1) - t73 * mrSges(4,3) + qJD(3) * t78 - t71 * t119 - t133;
t9 = m(4) * t122 - qJDD(3) * mrSges(4,2) + t74 * mrSges(4,3) - qJD(3) * t77 - t93 * t11 + t71 * t118 + t130 * t12;
t4 = m(3) * t114 - t99 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t94 * t10 + t96 * t9;
t6 = m(3) * t59 + t96 * t10 + t94 * t9;
t113 = m(2) * t88 + t4 * t127 + t90 * t131 + t92 * t6;
t103 = -t66 * mrSges(7,2) - t82 * t55 + t115;
t2 = m(2) * t76 + t97 * t4 - t95 * t8;
t1 = m(2) * t75 - t90 * t6 + (t4 * t95 + t131) * t92;
t3 = [-m(1) * g(1) - t89 * t1 + t91 * t2, t2, t4, t9, t12, -t42 * mrSges(6,2) - t69 * t47 + t103 + t137, t103; -m(1) * g(2) + t91 * t1 + t89 * t2, t1, t8, t10, t11, -t43 * mrSges(6,3) - t70 * t57 + t104, -t43 * mrSges(7,3) - t70 * t48 + t116; -m(1) * g(3) + t113, t113, t6, t134, t133, t66 * mrSges(6,1) + t82 * t58 + (t47 - t48) * t70 + (mrSges(6,2) - mrSges(7,3)) * t43 + t106, t110;];
f_new  = t3;
