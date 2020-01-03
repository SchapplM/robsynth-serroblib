% Calculate vector of cutting forces with Newton-Euler
% S5RRPPR8
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRPPR8_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR8_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR8_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR8_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR8_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR8_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR8_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:38:10
% EndTime: 2019-12-31 19:38:13
% DurationCPUTime: 1.17s
% Computational Cost: add. (9697->161), mult. (21904->210), div. (0->0), fcn. (12835->8), ass. (0->78)
t88 = sin(pkin(8));
t89 = cos(pkin(8));
t92 = sin(qJ(2));
t95 = cos(qJ(2));
t56 = (-t88 * t92 - t89 * t95) * qJD(1);
t57 = (-t88 * t95 + t89 * t92) * qJD(1);
t91 = sin(qJ(5));
t94 = cos(qJ(5));
t37 = t56 * t91 + t57 * t94;
t114 = qJD(1) * qJD(2);
t111 = t95 * t114;
t67 = qJDD(1) * t92 + t111;
t110 = t92 * t114;
t68 = qJDD(1) * t95 - t110;
t43 = -t67 * t88 - t68 * t89;
t44 = t67 * t89 - t68 * t88;
t18 = -t37 * qJD(5) + t43 * t94 - t44 * t91;
t36 = t56 * t94 - t57 * t91;
t19 = t36 * qJD(5) + t43 * t91 + t44 * t94;
t83 = -qJD(2) + qJD(5);
t34 = -mrSges(6,2) * t83 + t36 * mrSges(6,3);
t35 = mrSges(6,1) * t83 - t37 * mrSges(6,3);
t47 = -qJD(2) * pkin(4) - pkin(7) * t57;
t55 = t56 ^ 2;
t93 = sin(qJ(1));
t96 = cos(qJ(1));
t117 = t93 * g(1) - t96 * g(2);
t98 = qJD(1) ^ 2;
t58 = -qJDD(1) * pkin(1) - t98 * pkin(6) - t117;
t106 = -t68 * pkin(2) + t58 + (-t111 - t67) * qJ(3);
t116 = qJD(1) * t92;
t123 = t95 ^ 2 * t98;
t125 = 2 * qJD(3);
t70 = -qJD(2) * pkin(3) - qJ(4) * t116;
t99 = -pkin(2) * t110 + pkin(3) * t68 - qJ(4) * t123 + qJDD(4) - t106 + (t125 + t70) * t116;
t108 = m(6) * (-t43 * pkin(4) - t55 * pkin(7) + t57 * t47 + t99) + t19 * mrSges(6,2) - t18 * mrSges(6,1) + t37 * t35 - t36 * t34;
t45 = qJD(2) * mrSges(5,2) + mrSges(5,3) * t56;
t46 = -qJD(2) * mrSges(5,1) - mrSges(5,3) * t57;
t102 = m(5) * t99 - t43 * mrSges(5,1) + t44 * mrSges(5,2) - t56 * t45 + t57 * t46 + t108;
t101 = m(4) * ((pkin(2) * qJD(2) - (2 * qJD(3))) * t116 + t106) - t68 * mrSges(4,1) - t102;
t115 = qJD(1) * t95;
t74 = mrSges(4,2) * t115 + qJD(2) * mrSges(4,3);
t118 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t115 + t74;
t71 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t116;
t72 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t116;
t128 = -(t118 * t95 - (t71 - t72) * t92) * qJD(1) + m(3) * t58 - t68 * mrSges(3,1) + (mrSges(3,2) - mrSges(4,3)) * t67 + t101;
t107 = -g(1) * t96 - g(2) * t93;
t59 = -pkin(1) * t98 + qJDD(1) * pkin(6) + t107;
t112 = -g(3) * t92 + t95 * t59;
t64 = (-pkin(2) * t95 - qJ(3) * t92) * qJD(1);
t97 = qJD(2) ^ 2;
t103 = -pkin(2) * t97 + qJDD(2) * qJ(3) + qJD(2) * t125 + t64 * t115 + t112;
t65 = (-mrSges(4,1) * t95 - mrSges(4,3) * t92) * qJD(1);
t28 = -pkin(3) * t123 - qJ(4) * t68 + qJD(2) * t70 + t103;
t120 = -t95 * g(3) - t92 * t59;
t33 = -qJDD(2) * pkin(2) - qJ(3) * t97 + t64 * t116 + qJDD(3) - t120;
t29 = (-t67 + t111) * qJ(4) + (-t92 * t95 * t98 - qJDD(2)) * pkin(3) + t33;
t109 = -0.2e1 * qJD(4) * t57 - t88 * t28 + t89 * t29;
t12 = (-qJD(2) * t56 - t44) * pkin(7) + (t56 * t57 - qJDD(2)) * pkin(4) + t109;
t113 = 0.2e1 * qJD(4) * t56 + t89 * t28 + t88 * t29;
t13 = -pkin(4) * t55 + pkin(7) * t43 + qJD(2) * t47 + t113;
t27 = -mrSges(6,1) * t36 + mrSges(6,2) * t37;
t82 = -qJDD(2) + qJDD(5);
t10 = m(6) * (t12 * t94 - t13 * t91) - t19 * mrSges(6,3) + t82 * mrSges(6,1) - t37 * t27 + t83 * t34;
t11 = m(6) * (t12 * t91 + t13 * t94) + t18 * mrSges(6,3) - t82 * mrSges(6,2) + t36 * t27 - t83 * t35;
t40 = -mrSges(5,1) * t56 + mrSges(5,2) * t57;
t7 = m(5) * t109 - qJDD(2) * mrSges(5,1) - t44 * mrSges(5,3) - qJD(2) * t45 + t94 * t10 + t91 * t11 - t57 * t40;
t8 = m(5) * t113 + qJDD(2) * mrSges(5,2) + t43 * mrSges(5,3) + qJD(2) * t46 - t91 * t10 + t94 * t11 + t56 * t40;
t105 = m(4) * t103 + qJDD(2) * mrSges(4,3) + qJD(2) * t72 + t65 * t115 - t88 * t7 + t89 * t8;
t121 = mrSges(3,3) + mrSges(4,2);
t66 = (-mrSges(3,1) * t95 + mrSges(3,2) * t92) * qJD(1);
t4 = m(3) * t112 - qJDD(2) * mrSges(3,2) - qJD(2) * t71 + t66 * t115 + t121 * t68 + t105;
t104 = -m(4) * t33 - t89 * t7 - t88 * t8;
t5 = m(3) * t120 - t121 * t67 + (mrSges(3,1) + mrSges(4,1)) * qJDD(2) + t118 * qJD(2) + (-t65 - t66) * t116 + t104;
t124 = t92 * t4 + t95 * t5;
t9 = m(2) * t117 + qJDD(1) * mrSges(2,1) - t98 * mrSges(2,2) - t128;
t1 = m(2) * t107 - t98 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t95 * t4 - t92 * t5;
t2 = [-m(1) * g(1) + t1 * t96 - t9 * t93, t1, t4, t68 * mrSges(4,2) + t105, t8, t11; -m(1) * g(2) + t1 * t93 + t9 * t96, t9, t5, t101 + (-t92 * t72 - t95 * t74) * qJD(1) - t67 * mrSges(4,3), t7, t10; (-m(1) - m(2)) * g(3) + t124, -m(2) * g(3) + t124, t128, -qJDD(2) * mrSges(4,1) + t67 * mrSges(4,2) - qJD(2) * t74 + t65 * t116 - t104, t102, t108;];
f_new = t2;
