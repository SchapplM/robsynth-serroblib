% Calculate vector of cutting forces with Newton-Euler
% S5PRRPP4
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
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
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
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRRPP4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP4_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP4_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP4_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP4_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP4_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:40:48
% EndTime: 2019-12-31 17:40:49
% DurationCPUTime: 0.48s
% Computational Cost: add. (2148->133), mult. (4245->159), div. (0->0), fcn. (1994->6), ass. (0->63)
t64 = sin(qJ(3));
t66 = cos(qJ(3));
t83 = qJD(2) * qJD(3);
t39 = qJDD(2) * t64 + t66 * t83;
t40 = qJDD(2) * t66 - t64 * t83;
t85 = qJD(2) * t64;
t47 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t85;
t48 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t85;
t62 = sin(pkin(7));
t63 = cos(pkin(7));
t43 = g(1) * t62 - g(2) * t63;
t44 = -g(1) * t63 - g(2) * t62;
t65 = sin(qJ(2));
t67 = cos(qJ(2));
t91 = t67 * t43 - t65 * t44;
t76 = qJDD(2) * pkin(2) + t91;
t72 = -qJ(4) * t39 - t76;
t100 = 2 * qJD(4);
t45 = -qJD(3) * pkin(4) - qJ(5) * t85;
t46 = -qJD(3) * mrSges(6,1) - mrSges(6,3) * t85;
t84 = qJD(2) * t66;
t49 = qJD(3) * mrSges(6,2) - mrSges(6,3) * t84;
t60 = t66 ^ 2;
t69 = qJD(2) ^ 2;
t86 = qJ(4) * t66;
t99 = pkin(3) + pkin(4);
t77 = m(6) * (qJDD(5) + (-qJ(5) * t60 + pkin(6)) * t69 + t99 * t40 + (qJD(3) * t86 + (-pkin(3) * qJD(3) + t100 + t45) * t64) * qJD(2) - t72) + t39 * mrSges(6,2) + t40 * mrSges(6,1) + t46 * t85 + t49 * t84;
t97 = pkin(6) * t69;
t74 = m(5) * (-pkin(3) * t40 - t97 + (-0.2e1 * qJD(4) * t64 + (pkin(3) * t64 - t86) * qJD(3)) * qJD(2) + t72) - t40 * mrSges(5,1) - t77;
t51 = mrSges(5,2) * t84 + qJD(3) * mrSges(5,3);
t87 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t84 + t51;
t104 = (-t87 * t66 + (t47 - t48) * t64) * qJD(2) + (mrSges(4,2) - mrSges(5,3)) * t39 + m(4) * (-t76 - t97) - t40 * mrSges(4,1) + t74;
t36 = (-mrSges(5,1) * t66 - mrSges(5,3) * t64) * qJD(2);
t35 = (-pkin(3) * t66 - qJ(4) * t64) * qJD(2);
t68 = qJD(3) ^ 2;
t90 = t65 * t43 + t67 * t44;
t21 = -pkin(2) * t69 + qJDD(2) * pkin(6) + t90;
t61 = -g(3) + qJDD(1);
t92 = t66 * t21 + t64 * t61;
t73 = -pkin(3) * t68 + qJDD(3) * qJ(4) + qJD(3) * t100 + t35 * t84 + t92;
t102 = m(5) * t73 + qJDD(3) * mrSges(5,3) + qJD(3) * t48 + t36 * t84;
t98 = pkin(4) * t69;
t96 = t61 * t66;
t95 = -mrSges(5,1) - mrSges(6,1);
t93 = mrSges(4,3) + mrSges(5,2);
t37 = (mrSges(6,1) * t66 + mrSges(6,2) * t64) * qJD(2);
t89 = -t37 + (-mrSges(4,1) * t66 + mrSges(4,2) * t64) * qJD(2);
t81 = -0.2e1 * qJD(2) * qJD(5);
t78 = m(6) * (-qJ(5) * t40 + qJD(3) * t45 - t60 * t98 + t66 * t81 + t73) + qJD(3) * t46 + qJDD(3) * mrSges(6,2) - t40 * mrSges(6,3);
t7 = m(4) * t92 - qJDD(3) * mrSges(4,2) - qJD(3) * t47 + t93 * t40 + t89 * t84 + t102 + t78;
t18 = t64 * t21;
t75 = -qJ(4) * t68 + t35 * t85 + qJDD(4) + t18;
t12 = m(6) * (t64 * t81 - qJ(5) * t39 - t99 * qJDD(3) + (qJ(5) * t83 - t64 * t98 - t61) * t66 + t75);
t80 = m(5) * (-qJDD(3) * pkin(3) + t75 - t96) + t12;
t8 = m(4) * (-t18 + t96) + (mrSges(6,3) - t93) * t39 + (mrSges(4,1) - t95) * qJDD(3) + (t49 + t87) * qJD(3) + (-t36 - t89) * t85 - t80;
t82 = m(3) * t61 + t64 * t7 + t66 * t8;
t79 = m(2) * t61 + t82;
t71 = -t37 * t84 + t78;
t4 = m(3) * t91 + qJDD(2) * mrSges(3,1) - t69 * mrSges(3,2) - t104;
t3 = m(3) * t90 - t69 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t64 * t8 + t66 * t7;
t2 = m(2) * t44 + t3 * t67 - t4 * t65;
t1 = m(2) * t43 + t3 * t65 + t4 * t67;
t5 = [-m(1) * g(1) - t1 * t62 + t2 * t63, t2, t3, t7, t40 * mrSges(5,2) + t102 + t71, t71; -m(1) * g(2) + t1 * t63 + t2 * t62, t1, t4, t8, -t39 * mrSges(5,3) + (-t64 * t48 - t66 * t51) * qJD(2) + t74, -qJDD(3) * mrSges(6,1) - t39 * mrSges(6,3) - qJD(3) * t49 - t37 * t85 + t12; -m(1) * g(3) + t79, t79, t82, t104, (mrSges(5,2) - mrSges(6,3)) * t39 + t95 * qJDD(3) + (-t49 - t51) * qJD(3) + (t36 - t37) * t85 + t80, t77;];
f_new = t5;
