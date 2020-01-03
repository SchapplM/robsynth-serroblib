% Calculate vector of cutting forces with Newton-Euler
% S5RRPPR10
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRPPR10_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR10_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR10_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR10_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR10_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR10_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR10_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:43:05
% EndTime: 2019-12-31 19:43:07
% DurationCPUTime: 0.95s
% Computational Cost: add. (8510->163), mult. (18645->204), div. (0->0), fcn. (11421->8), ass. (0->80)
t83 = cos(qJ(2));
t105 = t83 * qJD(1);
t80 = sin(qJ(2));
t106 = qJD(1) * t80;
t107 = cos(pkin(8));
t78 = sin(pkin(8));
t60 = -t107 * qJD(2) + t78 * t106;
t102 = t60 * t105;
t104 = qJD(1) * qJD(2);
t99 = t83 * t104;
t66 = t80 * qJDD(1) + t99;
t49 = t78 * qJDD(2) + t107 * t66;
t86 = qJD(1) ^ 2;
t81 = sin(qJ(1));
t84 = cos(qJ(1));
t97 = -t84 * g(1) - t81 * g(2);
t57 = -t86 * pkin(1) + qJDD(1) * pkin(6) + t97;
t108 = -t83 * g(3) - t80 * t57;
t64 = (-pkin(2) * t83 - qJ(3) * t80) * qJD(1);
t85 = qJD(2) ^ 2;
t89 = qJDD(2) * pkin(2) + t85 * qJ(3) - t64 * t106 - qJDD(3) + t108;
t119 = -(t49 + t102) * qJ(4) - t89;
t116 = -2 * qJD(3);
t101 = t81 * g(1) - t84 * g(2);
t56 = -qJDD(1) * pkin(1) - t86 * pkin(6) - t101;
t72 = t80 * t104;
t67 = t83 * qJDD(1) - t72;
t28 = (-t66 - t99) * qJ(3) + (-t67 + t72) * pkin(2) + t56;
t100 = -t80 * g(3) + t83 * t57;
t32 = -t85 * pkin(2) + qJDD(2) * qJ(3) + t64 * t105 + t100;
t103 = t107 * t32 + t60 * t116 + t78 * t28;
t61 = t78 * qJD(2) + t107 * t106;
t47 = mrSges(5,1) * t105 + t61 * mrSges(5,2);
t109 = -mrSges(4,1) * t105 - t61 * mrSges(4,3) - t47;
t41 = t60 * mrSges(5,1) - t61 * mrSges(5,3);
t110 = -t60 * mrSges(4,1) - t61 * mrSges(4,2) - t41;
t111 = -mrSges(4,3) - mrSges(5,2);
t48 = -t107 * qJDD(2) + t78 * t66;
t113 = t83 ^ 2 * t86;
t40 = t60 * pkin(3) - t61 * qJ(4);
t94 = t107 * t28 - t78 * t32;
t17 = t67 * pkin(3) - qJ(4) * t113 + qJDD(4) - t94 + ((2 * qJD(3)) + t40) * t61;
t12 = (-t49 + t102) * pkin(7) + (t60 * t61 + t67) * pkin(4) + t17;
t50 = pkin(4) * t105 - t61 * pkin(7);
t58 = t60 ^ 2;
t115 = -2 * qJD(4);
t90 = -pkin(3) * t113 - t67 * qJ(4) + t105 * t115 - t60 * t40 + t103;
t13 = -t58 * pkin(4) + t48 * pkin(7) - t50 * t105 + t90;
t79 = sin(qJ(5));
t82 = cos(qJ(5));
t38 = t82 * t60 - t79 * t61;
t23 = t38 * qJD(5) + t79 * t48 + t82 * t49;
t39 = t79 * t60 + t82 * t61;
t26 = -t38 * mrSges(6,1) + t39 * mrSges(6,2);
t70 = qJD(5) + t105;
t33 = -t70 * mrSges(6,2) + t38 * mrSges(6,3);
t63 = qJDD(5) + t67;
t10 = m(6) * (t82 * t12 - t79 * t13) - t23 * mrSges(6,3) + t63 * mrSges(6,1) - t39 * t26 + t70 * t33;
t22 = -t39 * qJD(5) + t82 * t48 - t79 * t49;
t34 = t70 * mrSges(6,1) - t39 * mrSges(6,3);
t11 = m(6) * (t79 * t12 + t82 * t13) + t22 * mrSges(6,3) - t63 * mrSges(6,2) + t38 * t26 - t70 * t34;
t93 = m(5) * t90 - t67 * mrSges(5,3) - t79 * t10 + t82 * t11;
t5 = m(4) * t103 + t67 * mrSges(4,2) + t109 * t105 + t110 * t60 + t111 * t48 + t93;
t44 = -t60 * mrSges(5,2) - mrSges(5,3) * t105;
t45 = mrSges(4,2) * t105 - t60 * mrSges(4,3);
t92 = -m(5) * t17 - t82 * t10 - t79 * t11;
t6 = m(4) * t94 + (-mrSges(4,1) - mrSges(5,1)) * t67 + (m(4) * t116 + t110) * t61 + t111 * t49 + (-t44 - t45) * t105 + t92;
t68 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t106;
t69 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t105;
t118 = m(3) * t56 - t67 * mrSges(3,1) + t66 * mrSges(3,2) + (t68 * t80 - t69 * t83) * qJD(1) + t107 * t6 + t78 * t5;
t98 = m(6) * (-t58 * pkin(7) + (-pkin(3) - pkin(4)) * t48 + (pkin(3) * t105 + (2 * qJD(4)) + t50) * t61 - t119) + t23 * mrSges(6,2) - t22 * mrSges(6,1) + t39 * t34 - t38 * t33;
t91 = m(5) * (t61 * t115 + (-t61 * t105 + t48) * pkin(3) + t119) + t60 * t44 + t48 * mrSges(5,1) - t98;
t117 = -m(4) * t89 + t48 * mrSges(4,1) + t109 * t61 + (mrSges(4,2) - mrSges(5,3)) * t49 + t60 * t45 + t91;
t65 = (-mrSges(3,1) * t83 + mrSges(3,2) * t80) * qJD(1);
t4 = m(3) * t100 - qJDD(2) * mrSges(3,2) + t67 * mrSges(3,3) - qJD(2) * t68 + t65 * t105 + t107 * t5 - t78 * t6;
t8 = m(3) * t108 + qJDD(2) * mrSges(3,1) - t66 * mrSges(3,3) + qJD(2) * t69 - t65 * t106 - t117;
t114 = t80 * t4 + t83 * t8;
t2 = m(2) * t101 + qJDD(1) * mrSges(2,1) - t86 * mrSges(2,2) - t118;
t1 = m(2) * t97 - t86 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t83 * t4 - t80 * t8;
t3 = [-m(1) * g(1) + t84 * t1 - t81 * t2, t1, t4, t5, -t48 * mrSges(5,2) - t47 * t105 - t60 * t41 + t93, t11; -m(1) * g(2) + t81 * t1 + t84 * t2, t2, t8, t6, -t49 * mrSges(5,3) - t61 * t47 + t91, t10; (-m(1) - m(2)) * g(3) + t114, -m(2) * g(3) + t114, t118, t117, t67 * mrSges(5,1) + t49 * mrSges(5,2) + t44 * t105 + t61 * t41 - t92, t98;];
f_new = t3;
