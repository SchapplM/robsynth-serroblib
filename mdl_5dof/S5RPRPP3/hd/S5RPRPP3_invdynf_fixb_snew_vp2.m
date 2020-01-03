% Calculate vector of cutting forces with Newton-Euler
% S5RPRPP3
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
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
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRPP3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP3_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP3_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP3_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP3_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP3_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:22
% EndTime: 2019-12-31 18:12:23
% DurationCPUTime: 0.68s
% Computational Cost: add. (4068->151), mult. (9742->172), div. (0->0), fcn. (6060->6), ass. (0->70)
t72 = qJD(1) ^ 2;
t107 = cos(qJ(3));
t66 = sin(pkin(7));
t67 = cos(pkin(7));
t68 = sin(qJ(3));
t112 = -t67 * t107 + t66 * t68;
t80 = t107 * t66 + t67 * t68;
t55 = t80 * qJD(1);
t95 = t55 * qJD(3);
t39 = t112 * qJDD(1) + t95;
t54 = t112 * qJD(1);
t96 = t54 * qJD(3);
t40 = t80 * qJDD(1) - t96;
t42 = -qJD(3) * mrSges(4,2) - t54 * mrSges(4,3);
t43 = qJD(3) * mrSges(4,1) - t55 * mrSges(4,3);
t46 = t54 * mrSges(5,1) - qJD(3) * mrSges(5,3);
t69 = sin(qJ(1));
t70 = cos(qJ(1));
t90 = t69 * g(1) - t70 * g(2);
t86 = qJDD(2) - t90;
t64 = t67 ^ 2;
t98 = t66 ^ 2 + t64;
t75 = (-pkin(2) * t67 - pkin(1)) * qJDD(1) + (-t98 * pkin(6) - qJ(2)) * t72 + t86;
t48 = t55 * mrSges(5,1) + qJD(3) * mrSges(5,2);
t111 = -2 * qJD(4);
t73 = pkin(3) * t95 + t55 * t111 + (-t40 + t96) * qJ(4) + t75;
t44 = t55 * pkin(4) - qJD(3) * qJ(5);
t45 = t55 * mrSges(6,1) - qJD(3) * mrSges(6,3);
t47 = -t54 * mrSges(6,1) + qJD(3) * mrSges(6,2);
t53 = t54 ^ 2;
t79 = t40 * mrSges(6,2) + t55 * t45 - t54 * t47 - t39 * mrSges(6,3) - m(6) * (-t53 * pkin(4) + 0.2e1 * qJD(5) * t54 - t55 * t44 + (pkin(3) + qJ(5)) * t39 + t73);
t77 = -m(5) * (t39 * pkin(3) + t73) + t55 * t48 + t40 * mrSges(5,3) + t79;
t74 = m(4) * t75 + t40 * mrSges(4,2) + (t42 - t46) * t54 + (mrSges(4,1) - mrSges(5,2)) * t39 + t55 * t43 - t77;
t114 = -m(3) * (-qJDD(1) * pkin(1) - t72 * qJ(2) + t86) - t74;
t113 = t98 * mrSges(3,3);
t87 = -t70 * g(1) - t69 * g(2);
t56 = -t72 * pkin(1) + qJDD(1) * qJ(2) + t87;
t33 = -t54 * mrSges(5,2) - t55 * mrSges(5,3);
t101 = -t54 * mrSges(4,1) - t55 * mrSges(4,2) - t33;
t103 = -mrSges(4,3) - mrSges(5,1);
t104 = mrSges(5,2) - mrSges(6,3);
t31 = t54 * pkin(3) - t55 * qJ(4);
t71 = qJD(3) ^ 2;
t108 = pkin(2) * t72;
t94 = qJD(1) * qJD(2);
t89 = -t67 * g(3) - 0.2e1 * t66 * t94;
t97 = pkin(6) * qJDD(1);
t24 = (t67 * t108 - t56 - t97) * t66 + t89;
t88 = -t66 * g(3) + (t56 + 0.2e1 * t94) * t67;
t25 = -t64 * t108 + t67 * t97 + t88;
t85 = t107 * t24 - t68 * t25;
t18 = -qJDD(3) * pkin(3) - t71 * qJ(4) + t55 * t31 + qJDD(4) - t85;
t30 = -t55 * mrSges(6,2) + t54 * mrSges(6,3);
t93 = m(6) * (-0.2e1 * qJD(5) * qJD(3) + (t54 * t55 - qJDD(3)) * qJ(5) + (t40 + t96) * pkin(4) + t18) + t55 * t30 + t40 * mrSges(6,1);
t83 = m(5) * t18 + t93;
t99 = t46 - t47;
t7 = m(4) * t85 + t101 * t55 + t103 * t40 + (mrSges(4,1) - t104) * qJDD(3) + (t42 - t99) * qJD(3) - t83;
t102 = t107 * t25 + t68 * t24;
t78 = -t71 * pkin(3) + qJDD(3) * qJ(4) - t54 * t31 + t102;
t92 = m(6) * (-t39 * pkin(4) - t53 * qJ(5) + qJDD(5) + ((2 * qJD(4)) + t44) * qJD(3) + t78) + qJD(3) * t45 + qJDD(3) * mrSges(6,2);
t82 = m(5) * (qJD(3) * t111 - t78) - t92;
t8 = m(4) * t102 + (-mrSges(4,2) + mrSges(5,3)) * qJDD(3) + (-t43 + t48) * qJD(3) + (-t30 + t101) * t54 + (-mrSges(6,1) + t103) * t39 - t82;
t84 = -t67 * mrSges(3,1) + t66 * mrSges(3,2);
t81 = qJDD(1) * mrSges(3,3) + t72 * t84;
t4 = m(3) * t89 + t68 * t8 + t107 * t7 + (-m(3) * t56 - t81) * t66;
t5 = m(3) * t88 + t107 * t8 + t81 * t67 - t68 * t7;
t110 = t67 * t4 + t66 * t5;
t6 = (-mrSges(2,2) + t113) * t72 + (mrSges(2,1) - t84) * qJDD(1) + m(2) * t90 + t114;
t1 = m(2) * t87 - t72 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t66 * t4 + t67 * t5;
t2 = [-m(1) * g(1) + t70 * t1 - t69 * t6, t1, t5, t8, -t39 * mrSges(5,2) - t54 * t46 - t77, -t79; -m(1) * g(2) + t69 * t1 + t70 * t6, t6, t4, t7, -qJDD(3) * mrSges(5,3) - qJD(3) * t48 + (t30 + t33) * t54 + (mrSges(5,1) + mrSges(6,1)) * t39 + t82, -qJDD(3) * mrSges(6,3) - qJD(3) * t47 + t93; (-m(1) - m(2)) * g(3) + t110, -m(2) * g(3) + t110, t84 * qJDD(1) - t72 * t113 - t114, t74, t40 * mrSges(5,1) + t99 * qJD(3) + t104 * qJDD(3) + t55 * t33 + t83, -t39 * mrSges(6,1) - t54 * t30 + t92;];
f_new = t2;
