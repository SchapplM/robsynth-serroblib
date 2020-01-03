% Calculate vector of cutting forces with Newton-Euler
% S5RRPRR8
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRPRR8_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR8_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR8_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR8_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR8_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR8_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR8_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:16:59
% EndTime: 2019-12-31 20:17:03
% DurationCPUTime: 2.12s
% Computational Cost: add. (24341->166), mult. (56452->224), div. (0->0), fcn. (40019->10), ass. (0->85)
t104 = qJD(1) * qJD(2);
t83 = sin(qJ(2));
t87 = cos(qJ(2));
t69 = qJDD(1) * t83 + t104 * t87;
t70 = qJDD(1) * t87 - t104 * t83;
t106 = qJD(1) * t83;
t72 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t106;
t105 = qJD(1) * t87;
t73 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t105;
t89 = qJD(1) ^ 2;
t79 = sin(pkin(9));
t80 = cos(pkin(9));
t54 = -t69 * t79 + t70 * t80;
t55 = t69 * t80 + t70 * t79;
t63 = (-t79 * t83 + t80 * t87) * qJD(1);
t56 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t63;
t64 = (t79 * t87 + t80 * t83) * qJD(1);
t57 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t64;
t71 = qJD(2) * pkin(2) - qJ(3) * t106;
t78 = t87 ^ 2;
t84 = sin(qJ(1));
t88 = cos(qJ(1));
t102 = t84 * g(1) - t88 * g(2);
t96 = -qJDD(1) * pkin(1) - t102;
t93 = -t70 * pkin(2) + qJDD(3) + t71 * t106 + (-qJ(3) * t78 - pkin(6)) * t89 + t96;
t99 = -g(1) * t88 - g(2) * t84;
t66 = -pkin(1) * t89 + qJDD(1) * pkin(6) + t99;
t108 = t83 * t66;
t109 = pkin(2) * t89;
t39 = qJDD(2) * pkin(2) - t69 * qJ(3) - t108 + (qJ(3) * t104 + t109 * t83 - g(3)) * t87;
t101 = -g(3) * t83 + t87 * t66;
t40 = t70 * qJ(3) - qJD(2) * t71 - t109 * t78 + t101;
t100 = -0.2e1 * qJD(3) * t64 + t80 * t39 - t79 * t40;
t19 = (qJD(2) * t63 - t55) * pkin(7) + (t63 * t64 + qJDD(2)) * pkin(3) + t100;
t103 = 0.2e1 * qJD(3) * t63 + t79 * t39 + t80 * t40;
t58 = qJD(2) * pkin(3) - pkin(7) * t64;
t62 = t63 ^ 2;
t21 = -pkin(3) * t62 + pkin(7) * t54 - qJD(2) * t58 + t103;
t82 = sin(qJ(4));
t86 = cos(qJ(4));
t107 = t19 * t82 + t86 * t21;
t48 = t63 * t86 - t64 * t82;
t49 = t63 * t82 + t64 * t86;
t35 = -pkin(4) * t48 - pkin(8) * t49;
t77 = qJD(2) + qJD(4);
t75 = t77 ^ 2;
t76 = qJDD(2) + qJDD(4);
t16 = -pkin(4) * t75 + pkin(8) * t76 + t48 * t35 + t107;
t28 = -t49 * qJD(4) + t54 * t86 - t55 * t82;
t29 = t48 * qJD(4) + t54 * t82 + t55 * t86;
t91 = -t54 * pkin(3) - t62 * pkin(7) + t64 * t58 + t93;
t17 = (-t48 * t77 - t29) * pkin(8) + (t49 * t77 - t28) * pkin(4) + t91;
t81 = sin(qJ(5));
t85 = cos(qJ(5));
t43 = -t49 * t81 + t77 * t85;
t23 = t43 * qJD(5) + t29 * t85 + t76 * t81;
t27 = qJDD(5) - t28;
t44 = t49 * t85 + t77 * t81;
t30 = -mrSges(6,1) * t43 + mrSges(6,2) * t44;
t47 = qJD(5) - t48;
t31 = -mrSges(6,2) * t47 + mrSges(6,3) * t43;
t13 = m(6) * (-t16 * t81 + t17 * t85) - t23 * mrSges(6,3) + t27 * mrSges(6,1) - t44 * t30 + t47 * t31;
t22 = -t44 * qJD(5) - t29 * t81 + t76 * t85;
t32 = mrSges(6,1) * t47 - mrSges(6,3) * t44;
t14 = m(6) * (t16 * t85 + t17 * t81) + t22 * mrSges(6,3) - t27 * mrSges(6,2) + t43 * t30 - t47 * t32;
t45 = -mrSges(5,2) * t77 + t48 * mrSges(5,3);
t46 = mrSges(5,1) * t77 - t49 * mrSges(5,3);
t95 = -m(5) * t91 + t28 * mrSges(5,1) - t29 * mrSges(5,2) - t13 * t85 - t14 * t81 + t48 * t45 - t49 * t46;
t92 = -m(4) * t93 + t54 * mrSges(4,1) - t55 * mrSges(4,2) + t63 * t56 - t64 * t57 + t95;
t111 = (t72 * t83 - t73 * t87) * qJD(1) + m(3) * (-t89 * pkin(6) + t96) - t70 * mrSges(3,1) + t69 * mrSges(3,2) - t92;
t34 = -mrSges(5,1) * t48 + mrSges(5,2) * t49;
t98 = t19 * t86 - t21 * t82;
t94 = m(6) * (-pkin(4) * t76 - pkin(8) * t75 + t49 * t35 - t98) - t22 * mrSges(6,1) + t23 * mrSges(6,2) - t43 * t31 + t44 * t32;
t10 = m(5) * t98 + mrSges(5,1) * t76 - t29 * mrSges(5,3) - t49 * t34 + t45 * t77 - t94;
t52 = -mrSges(4,1) * t63 + mrSges(4,2) * t64;
t9 = m(5) * t107 - mrSges(5,2) * t76 + t28 * mrSges(5,3) - t13 * t81 + t14 * t85 + t48 * t34 - t46 * t77;
t6 = m(4) * t100 + qJDD(2) * mrSges(4,1) - t55 * mrSges(4,3) + qJD(2) * t56 + t10 * t86 - t64 * t52 + t82 * t9;
t68 = (-mrSges(3,1) * t87 + mrSges(3,2) * t83) * qJD(1);
t7 = m(4) * t103 - qJDD(2) * mrSges(4,2) + t54 * mrSges(4,3) - qJD(2) * t57 - t10 * t82 + t63 * t52 + t86 * t9;
t4 = m(3) * (-t87 * g(3) - t108) - t69 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t68 * t106 + qJD(2) * t73 + t79 * t7 + t80 * t6;
t5 = m(3) * t101 - qJDD(2) * mrSges(3,2) + t70 * mrSges(3,3) - qJD(2) * t72 + t105 * t68 - t6 * t79 + t7 * t80;
t110 = t4 * t87 + t5 * t83;
t8 = m(2) * t102 + qJDD(1) * mrSges(2,1) - t89 * mrSges(2,2) - t111;
t1 = m(2) * t99 - mrSges(2,1) * t89 - qJDD(1) * mrSges(2,2) - t4 * t83 + t5 * t87;
t2 = [-m(1) * g(1) + t1 * t88 - t8 * t84, t1, t5, t7, t9, t14; -m(1) * g(2) + t1 * t84 + t8 * t88, t8, t4, t6, t10, t13; (-m(1) - m(2)) * g(3) + t110, -m(2) * g(3) + t110, t111, -t92, -t95, t94;];
f_new = t2;
