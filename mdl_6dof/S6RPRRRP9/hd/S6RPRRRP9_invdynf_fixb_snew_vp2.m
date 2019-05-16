% Calculate vector of cutting forces with Newton-Euler
% S6RPRRRP9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2019-05-06 01:54
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRRRP9_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP9_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP9_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP9_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP9_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP9_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP9_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:50:49
% EndTime: 2019-05-06 01:50:53
% DurationCPUTime: 1.56s
% Computational Cost: add. (18558->175), mult. (36185->212), div. (0->0), fcn. (23148->8), ass. (0->86)
t91 = cos(qJ(3));
t119 = qJD(1) * t91;
t86 = sin(qJ(4));
t90 = cos(qJ(4));
t70 = t86 * qJD(3) + t90 * t119;
t118 = qJD(1) * qJD(3);
t87 = sin(qJ(3));
t108 = t87 * t118;
t74 = t91 * qJDD(1) - t108;
t46 = -t70 * qJD(4) + t90 * qJDD(3) - t86 * t74;
t69 = t90 * qJD(3) - t86 * t119;
t47 = t69 * qJD(4) + t86 * qJDD(3) + t90 * t74;
t85 = sin(qJ(5));
t89 = cos(qJ(5));
t50 = t85 * t69 + t89 * t70;
t26 = -t50 * qJD(5) + t89 * t46 - t85 * t47;
t49 = t89 * t69 - t85 * t70;
t27 = t49 * qJD(5) + t85 * t46 + t89 * t47;
t81 = t87 * qJD(1);
t78 = t81 + qJD(4);
t77 = qJD(5) + t78;
t42 = t77 * pkin(5) - t50 * qJ(6);
t43 = t77 * mrSges(7,1) - t50 * mrSges(7,3);
t48 = t49 ^ 2;
t88 = sin(qJ(1));
t92 = cos(qJ(1));
t113 = t88 * g(1) - t92 * g(2);
t94 = qJD(1) ^ 2;
t100 = -t94 * qJ(2) + qJDD(2) - t113;
t126 = -pkin(1) - pkin(7);
t58 = t126 * qJDD(1) + t100;
t104 = t87 * g(3) + t91 * t58;
t72 = (pkin(3) * t87 - pkin(8) * t91) * qJD(1);
t93 = qJD(3) ^ 2;
t38 = -qJDD(3) * pkin(3) - t93 * pkin(8) + t72 * t119 - t104;
t56 = t78 * pkin(4) - t70 * pkin(9);
t67 = t69 ^ 2;
t96 = -t46 * pkin(4) - t67 * pkin(9) + t70 * t56 + t38;
t114 = m(7) * (-t26 * pkin(5) - t48 * qJ(6) + t50 * t42 + qJDD(6) + t96) + t27 * mrSges(7,2) + t50 * t43;
t40 = -t77 * mrSges(7,2) + t49 * mrSges(7,3);
t41 = -t77 * mrSges(6,2) + t49 * mrSges(6,3);
t44 = t77 * mrSges(6,1) - t50 * mrSges(6,3);
t131 = m(6) * t96 + t27 * mrSges(6,2) + t50 * t44 + t114 - (t41 + t40) * t49 - (mrSges(6,1) + mrSges(7,1)) * t26;
t52 = -t78 * mrSges(5,2) + t69 * mrSges(5,3);
t53 = t78 * mrSges(5,1) - t70 * mrSges(5,3);
t130 = m(5) * t38 - t46 * mrSges(5,1) + t47 * mrSges(5,2) - t69 * t52 + t70 * t53 + t131;
t105 = -t92 * g(1) - t88 * g(2);
t129 = -qJDD(1) * qJ(2) - 0.2e1 * qJD(2) * qJD(1) - t105;
t127 = -m(2) - m(3);
t125 = mrSges(2,1) - mrSges(3,2);
t123 = -mrSges(2,2) + mrSges(3,3);
t79 = t91 * t118;
t73 = -t87 * qJDD(1) - t79;
t97 = t126 * t94 - t129;
t35 = (-t74 + t108) * pkin(8) + (-t73 + t79) * pkin(3) + t97;
t112 = -t91 * g(3) + t87 * t58;
t39 = -t93 * pkin(3) + qJDD(3) * pkin(8) - t72 * t81 + t112;
t106 = t90 * t35 - t86 * t39;
t68 = qJDD(4) - t73;
t18 = (t69 * t78 - t47) * pkin(9) + (t69 * t70 + t68) * pkin(4) + t106;
t121 = t86 * t35 + t90 * t39;
t20 = -t67 * pkin(4) + t46 * pkin(9) - t78 * t56 + t121;
t122 = t85 * t18 + t89 * t20;
t107 = t89 * t18 - t85 * t20;
t63 = qJDD(5) + t68;
t116 = m(7) * (-0.2e1 * qJD(6) * t50 + (t49 * t77 - t27) * qJ(6) + (t49 * t50 + t63) * pkin(5) + t107) + t77 * t40 + t63 * mrSges(7,1);
t31 = -t49 * mrSges(7,1) + t50 * mrSges(7,2);
t115 = m(7) * (-t48 * pkin(5) + t26 * qJ(6) + 0.2e1 * qJD(6) * t49 - t77 * t42 + t122) + t26 * mrSges(7,3) + t49 * t31;
t71 = (mrSges(4,1) * t87 + mrSges(4,2) * t91) * qJD(1);
t75 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t81;
t11 = m(4) * t104 + qJDD(3) * mrSges(4,1) - t74 * mrSges(4,3) + qJD(3) * t75 - t71 * t119 - t130;
t32 = -t49 * mrSges(6,1) + t50 * mrSges(6,2);
t10 = m(6) * t122 + t26 * mrSges(6,3) + t49 * t32 + (-t44 - t43) * t77 + (-mrSges(6,2) - mrSges(7,2)) * t63 + t115;
t51 = -t69 * mrSges(5,1) + t70 * mrSges(5,2);
t9 = m(6) * t107 + t63 * mrSges(6,1) + t77 * t41 + (-t32 - t31) * t50 + (-mrSges(6,3) - mrSges(7,3)) * t27 + t116;
t7 = m(5) * t106 + t68 * mrSges(5,1) - t47 * mrSges(5,3) + t85 * t10 - t70 * t51 + t78 * t52 + t89 * t9;
t76 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t119;
t8 = m(5) * t121 - t68 * mrSges(5,2) + t46 * mrSges(5,3) + t89 * t10 + t69 * t51 - t78 * t53 - t85 * t9;
t4 = m(4) * t112 - qJDD(3) * mrSges(4,2) + t73 * mrSges(4,3) - qJD(3) * t76 - t86 * t7 - t71 * t81 + t90 * t8;
t111 = -t87 * t11 + t91 * t4;
t101 = -m(3) * (-qJDD(1) * pkin(1) + t100) - t91 * t11 - t87 * t4;
t99 = m(4) * t97 - t73 * mrSges(4,1) + t74 * mrSges(4,2) + t76 * t119 + t90 * t7 + t75 * t81 + t86 * t8;
t98 = -m(3) * (t94 * pkin(1) + t129) + t99;
t2 = m(2) * t105 + t123 * qJDD(1) - t125 * t94 + t98;
t1 = m(2) * t113 + t125 * qJDD(1) + t123 * t94 + t101;
t3 = [-m(1) * g(1) - t88 * t1 + t92 * t2, t2, -m(3) * g(3) + t111, t4, t8, t10, -t63 * mrSges(7,2) - t77 * t43 + t115; -m(1) * g(2) + t92 * t1 + t88 * t2, t1, -t94 * mrSges(3,2) - qJDD(1) * mrSges(3,3) - t98, t11, t7, t9, -t27 * mrSges(7,3) - t50 * t31 + t116; (-m(1) + t127) * g(3) + t111, t127 * g(3) + t111, qJDD(1) * mrSges(3,2) - t94 * mrSges(3,3) - t101, t99, t130, t131, -t26 * mrSges(7,1) - t49 * t40 + t114;];
f_new  = t3;
