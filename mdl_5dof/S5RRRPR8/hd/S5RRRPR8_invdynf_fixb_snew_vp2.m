% Calculate vector of cutting forces with Newton-Euler
% S5RRRPR8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRPR8_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR8_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR8_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR8_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR8_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR8_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR8_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:19:10
% EndTime: 2019-12-31 21:19:13
% DurationCPUTime: 1.08s
% Computational Cost: add. (9580->164), mult. (19854->206), div. (0->0), fcn. (12828->8), ass. (0->81)
t77 = sin(qJ(2));
t80 = cos(qJ(2));
t98 = qJD(1) * qJD(2);
t64 = t77 * qJDD(1) + t80 * t98;
t65 = t80 * qJDD(1) - t77 * t98;
t100 = qJD(1) * t77;
t66 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t100;
t99 = qJD(1) * t80;
t67 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t99;
t82 = qJD(1) ^ 2;
t110 = cos(qJ(3));
t76 = sin(qJ(3));
t57 = t76 * t100 - t110 * t99;
t73 = qJD(2) + qJD(3);
t49 = t57 * mrSges(5,1) - t73 * mrSges(5,3);
t101 = -t73 * mrSges(4,2) - t57 * mrSges(4,3) - t49;
t105 = mrSges(4,1) - mrSges(5,2);
t58 = (t110 * t77 + t76 * t80) * qJD(1);
t35 = t58 * qJD(3) - t110 * t65 + t76 * t64;
t36 = -t57 * qJD(3) + t110 * t64 + t76 * t65;
t48 = t73 * mrSges(4,1) - t58 * mrSges(4,3);
t68 = qJD(2) * pkin(2) - pkin(7) * t100;
t74 = t80 ^ 2;
t78 = sin(qJ(1));
t81 = cos(qJ(1));
t97 = t78 * g(1) - t81 * g(2);
t92 = -qJDD(1) * pkin(1) - t97;
t85 = -t65 * pkin(2) + t68 * t100 + (-pkin(7) * t74 - pkin(6)) * t82 + t92;
t51 = t58 * pkin(4) - t73 * pkin(8);
t56 = t57 ^ 2;
t109 = t57 * t73;
t114 = -2 * qJD(4);
t83 = (-t36 + t109) * qJ(4) + t85 + (t73 * pkin(3) + t114) * t58;
t12 = -t56 * pkin(4) - t58 * t51 + (pkin(3) + pkin(8)) * t35 + t83;
t42 = t57 * pkin(3) - t58 * qJ(4);
t71 = t73 ^ 2;
t72 = qJDD(2) + qJDD(3);
t95 = -t81 * g(1) - t78 * g(2);
t60 = -t82 * pkin(1) + qJDD(1) * pkin(6) + t95;
t106 = t77 * t60;
t111 = pkin(2) * t82;
t28 = qJDD(2) * pkin(2) - t64 * pkin(7) - t106 + (pkin(7) * t98 + t77 * t111 - g(3)) * t80;
t96 = -t77 * g(3) + t80 * t60;
t29 = t65 * pkin(7) - qJD(2) * t68 - t74 * t111 + t96;
t94 = t110 * t28 - t76 * t29;
t19 = -t72 * pkin(3) - t71 * qJ(4) + t58 * t42 + qJDD(4) - t94;
t13 = (t57 * t58 - t72) * pkin(8) + (t36 + t109) * pkin(4) + t19;
t75 = sin(qJ(5));
t79 = cos(qJ(5));
t45 = t79 * t57 - t75 * t73;
t22 = t45 * qJD(5) + t75 * t35 + t79 * t72;
t46 = t75 * t57 + t79 * t73;
t27 = -t45 * mrSges(6,1) + t46 * mrSges(6,2);
t33 = qJDD(5) + t36;
t55 = qJD(5) + t58;
t38 = -t55 * mrSges(6,2) + t45 * mrSges(6,3);
t10 = m(6) * (-t75 * t12 + t79 * t13) - t22 * mrSges(6,3) + t33 * mrSges(6,1) - t46 * t27 + t55 * t38;
t21 = -t46 * qJD(5) + t79 * t35 - t75 * t72;
t39 = t55 * mrSges(6,1) - t46 * mrSges(6,3);
t11 = m(6) * (t79 * t12 + t75 * t13) + t21 * mrSges(6,3) - t33 * mrSges(6,2) + t45 * t27 - t55 * t39;
t50 = t58 * mrSges(5,1) + t73 * mrSges(5,2);
t91 = t75 * t10 - t79 * t11 - m(5) * (t35 * pkin(3) + t83) + t36 * mrSges(5,3) + t58 * t50;
t84 = m(4) * t85 + t36 * mrSges(4,2) + t101 * t57 + t105 * t35 + t58 * t48 - t91;
t116 = (t77 * t66 - t80 * t67) * qJD(1) - t65 * mrSges(3,1) + t64 * mrSges(3,2) + m(3) * (-t82 * pkin(6) + t92) + t84;
t63 = (-mrSges(3,1) * t80 + mrSges(3,2) * t77) * qJD(1);
t44 = -t57 * mrSges(5,2) - t58 * mrSges(5,3);
t102 = -t57 * mrSges(4,1) - t58 * mrSges(4,2) - t44;
t104 = -mrSges(4,3) - mrSges(5,1);
t90 = -m(5) * t19 - t79 * t10 - t75 * t11;
t7 = m(4) * t94 + t101 * t73 + t102 * t58 + t104 * t36 + t105 * t72 + t90;
t103 = t110 * t29 + t76 * t28;
t87 = -t71 * pkin(3) + t72 * qJ(4) - t57 * t42 + t103;
t89 = -t21 * mrSges(6,1) - t45 * t38 + m(6) * (-t35 * pkin(4) - t56 * pkin(8) + ((2 * qJD(4)) + t51) * t73 + t87) + t22 * mrSges(6,2) + t46 * t39;
t86 = -m(5) * (t73 * t114 - t87) + t89;
t8 = m(4) * t103 + (-t48 + t50) * t73 + (-mrSges(4,2) + mrSges(5,3)) * t72 + t102 * t57 + t104 * t35 + t86;
t4 = m(3) * (-t80 * g(3) - t106) - t64 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t63 * t100 + qJD(2) * t67 + t76 * t8 + t110 * t7;
t5 = m(3) * t96 - qJDD(2) * mrSges(3,2) + t65 * mrSges(3,3) - qJD(2) * t66 + t110 * t8 + t63 * t99 - t76 * t7;
t113 = t80 * t4 + t77 * t5;
t6 = m(2) * t97 + qJDD(1) * mrSges(2,1) - t82 * mrSges(2,2) - t116;
t1 = m(2) * t95 - t82 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t77 * t4 + t80 * t5;
t2 = [-m(1) * g(1) + t81 * t1 - t78 * t6, t1, t5, t8, -t35 * mrSges(5,2) - t57 * t49 - t91, t11; -m(1) * g(2) + t78 * t1 + t81 * t6, t6, t4, t7, t35 * mrSges(5,1) - t72 * mrSges(5,3) + t57 * t44 - t73 * t50 - t86, t10; (-m(1) - m(2)) * g(3) + t113, -m(2) * g(3) + t113, t116, t84, t36 * mrSges(5,1) + t72 * mrSges(5,2) + t58 * t44 + t73 * t49 - t90, t89;];
f_new = t2;
