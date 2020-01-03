% Calculate vector of cutting forces with Newton-Euler
% S5RPRPP2
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
% Datum: 2019-12-31 18:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRPP2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP2_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP2_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP2_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP2_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP2_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP2_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP2_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:10:42
% EndTime: 2019-12-31 18:10:43
% DurationCPUTime: 0.50s
% Computational Cost: add. (2384->140), mult. (4602->165), div. (0->0), fcn. (1994->6), ass. (0->64)
t64 = sin(qJ(3));
t66 = cos(qJ(3));
t84 = qJD(1) * qJD(3);
t41 = qJDD(1) * t64 + t66 * t84;
t42 = qJDD(1) * t66 - t64 * t84;
t86 = qJD(1) * t64;
t47 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t86;
t48 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t86;
t65 = sin(qJ(1));
t67 = cos(qJ(1));
t80 = t65 * g(1) - g(2) * t67;
t35 = qJDD(1) * pkin(1) + t80;
t69 = qJD(1) ^ 2;
t77 = -g(1) * t67 - g(2) * t65;
t40 = -pkin(1) * t69 + t77;
t62 = sin(pkin(7));
t63 = cos(pkin(7));
t92 = t35 * t63 - t40 * t62;
t76 = qJDD(1) * pkin(2) + t92;
t72 = -qJ(4) * t41 - t76;
t100 = pkin(3) + pkin(4);
t101 = 2 * qJD(4);
t45 = -qJD(3) * pkin(4) - qJ(5) * t86;
t46 = -qJD(3) * mrSges(6,1) - mrSges(6,3) * t86;
t85 = qJD(1) * t66;
t49 = qJD(3) * mrSges(6,2) - mrSges(6,3) * t85;
t60 = t66 ^ 2;
t87 = qJ(4) * t66;
t78 = m(6) * (qJDD(5) + (-qJ(5) * t60 + pkin(6)) * t69 + t100 * t42 + (qJD(3) * t87 + (-pkin(3) * qJD(3) + t101 + t45) * t64) * qJD(1) - t72) + t41 * mrSges(6,2) + t42 * mrSges(6,1) + t46 * t86 + t49 * t85;
t98 = pkin(6) * t69;
t74 = m(5) * (-pkin(3) * t42 - t98 + (-0.2e1 * qJD(4) * t64 + (pkin(3) * t64 - t87) * qJD(3)) * qJD(1) + t72) - t42 * mrSges(5,1) - t78;
t51 = mrSges(5,2) * t85 + qJD(3) * mrSges(5,3);
t88 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t85 + t51;
t105 = (-t88 * t66 + (t47 - t48) * t64) * qJD(1) + (mrSges(4,2) - mrSges(5,3)) * t41 + m(4) * (-t76 - t98) - t42 * mrSges(4,1) + t74;
t37 = (-mrSges(5,1) * t66 - mrSges(5,3) * t64) * qJD(1);
t36 = (-pkin(3) * t66 - qJ(4) * t64) * qJD(1);
t68 = qJD(3) ^ 2;
t91 = t35 * t62 + t40 * t63;
t21 = -pkin(2) * t69 + qJDD(1) * pkin(6) + t91;
t61 = -g(3) + qJDD(2);
t93 = t21 * t66 + t64 * t61;
t73 = -pkin(3) * t68 + qJDD(3) * qJ(4) + qJD(3) * t101 + t36 * t85 + t93;
t103 = m(5) * t73 + qJDD(3) * mrSges(5,3) + qJD(3) * t48 + t37 * t85;
t99 = pkin(4) * t69;
t97 = t61 * t66;
t96 = -mrSges(5,1) - mrSges(6,1);
t94 = mrSges(4,3) + mrSges(5,2);
t38 = (mrSges(6,1) * t66 + mrSges(6,2) * t64) * qJD(1);
t90 = -t38 + (-mrSges(4,1) * t66 + mrSges(4,2) * t64) * qJD(1);
t82 = -0.2e1 * qJD(1) * qJD(5);
t79 = m(6) * (-qJ(5) * t42 + qJD(3) * t45 - t60 * t99 + t66 * t82 + t73) + qJD(3) * t46 + qJDD(3) * mrSges(6,2) - t42 * mrSges(6,3);
t7 = m(4) * t93 - qJDD(3) * mrSges(4,2) - qJD(3) * t47 + t42 * t94 + t85 * t90 + t103 + t79;
t18 = t64 * t21;
t75 = -qJ(4) * t68 + t36 * t86 + qJDD(4) + t18;
t12 = m(6) * (t64 * t82 - qJ(5) * t41 - t100 * qJDD(3) + (qJ(5) * t84 - t64 * t99 - t61) * t66 + t75);
t81 = m(5) * (-qJDD(3) * pkin(3) + t75 - t97) + t12;
t8 = m(4) * (-t18 + t97) + (mrSges(6,3) - t94) * t41 + (mrSges(4,1) - t96) * qJDD(3) + (t49 + t88) * qJD(3) + (-t37 - t90) * t86 - t81;
t83 = m(3) * t61 + t64 * t7 + t66 * t8;
t71 = -t38 * t85 + t79;
t4 = m(3) * t92 + qJDD(1) * mrSges(3,1) - t69 * mrSges(3,2) - t105;
t3 = m(3) * t91 - mrSges(3,1) * t69 - qJDD(1) * mrSges(3,2) - t64 * t8 + t66 * t7;
t2 = m(2) * t77 - mrSges(2,1) * t69 - qJDD(1) * mrSges(2,2) + t3 * t63 - t4 * t62;
t1 = m(2) * t80 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t69 + t3 * t62 + t4 * t63;
t5 = [-m(1) * g(1) - t1 * t65 + t2 * t67, t2, t3, t7, t42 * mrSges(5,2) + t103 + t71, t71; -m(1) * g(2) + t1 * t67 + t2 * t65, t1, t4, t8, -t41 * mrSges(5,3) + (-t64 * t48 - t66 * t51) * qJD(1) + t74, -qJDD(3) * mrSges(6,1) - t41 * mrSges(6,3) - qJD(3) * t49 - t38 * t86 + t12; (-m(1) - m(2)) * g(3) + t83, -m(2) * g(3) + t83, t83, t105, (mrSges(5,2) - mrSges(6,3)) * t41 + t96 * qJDD(3) + (-t49 - t51) * qJD(3) + (t37 - t38) * t86 + t81, t78;];
f_new = t5;
