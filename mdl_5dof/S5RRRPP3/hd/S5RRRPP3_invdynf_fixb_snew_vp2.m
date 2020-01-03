% Calculate vector of cutting forces with Newton-Euler
% S5RRRPP3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 20:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRPP3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP3_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP3_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP3_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP3_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP3_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:53:09
% EndTime: 2019-12-31 20:53:11
% DurationCPUTime: 0.53s
% Computational Cost: add. (3630->143), mult. (4580->163), div. (0->0), fcn. (1978->6), ass. (0->68)
t56 = qJD(1) + qJD(2);
t54 = t56 ^ 2;
t55 = qJDD(1) + qJDD(2);
t61 = sin(qJ(1));
t64 = cos(qJ(1));
t77 = t61 * g(1) - t64 * g(2);
t40 = qJDD(1) * pkin(1) + t77;
t66 = qJD(1) ^ 2;
t72 = -t64 * g(1) - t61 * g(2);
t41 = -t66 * pkin(1) + t72;
t60 = sin(qJ(2));
t63 = cos(qJ(2));
t75 = t63 * t40 - t60 * t41;
t70 = -t55 * pkin(2) - t75;
t21 = -t54 * pkin(7) + t70;
t59 = sin(qJ(3));
t62 = cos(qJ(3));
t82 = qJD(3) * t56;
t78 = t62 * t82;
t34 = t59 * t55 + t78;
t79 = t59 * t82;
t35 = t62 * t55 - t79;
t95 = t56 * t59;
t42 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t95;
t94 = t56 * t62;
t43 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t94;
t47 = mrSges(6,1) * t94 + qJD(3) * mrSges(6,2);
t101 = -2 * qJD(4);
t104 = pkin(3) * t79 + t95 * t101;
t46 = -mrSges(5,1) * t94 - qJD(3) * mrSges(5,3);
t100 = -2 * qJD(5);
t44 = pkin(4) * t95 - qJD(3) * qJ(5);
t58 = t62 ^ 2;
t97 = m(6) * (-t34 * qJ(4) + (-pkin(4) * t58 - pkin(7)) * t54 + (-pkin(3) - qJ(5)) * t35 + (-t44 * t59 + (-qJ(4) * qJD(3) + t100) * t62) * t56 + t70 + t104) - t35 * mrSges(6,3);
t73 = m(5) * (-t35 * pkin(3) + (-t34 - t78) * qJ(4) + t21 + t104) + t46 * t94 - t34 * mrSges(5,3) + t97;
t45 = mrSges(6,1) * t95 - qJD(3) * mrSges(6,3);
t48 = mrSges(5,1) * t95 + qJD(3) * mrSges(5,2);
t84 = -t45 - t48;
t105 = ((t42 + t84) * t59 - (t43 + t47) * t62) * t56 + m(4) * t21 + (mrSges(4,2) - mrSges(6,2)) * t34 - (mrSges(4,1) - mrSges(5,2)) * t35 + t73;
t99 = -m(2) - m(3);
t30 = (-pkin(3) * t62 - qJ(4) * t59) * t56;
t65 = qJD(3) ^ 2;
t86 = t60 * t40 + t63 * t41;
t22 = -t54 * pkin(2) + t55 * pkin(7) + t86;
t76 = -t59 * g(3) + t62 * t22;
t67 = -t65 * pkin(3) + qJDD(3) * qJ(4) + t30 * t94 + t76;
t33 = (-mrSges(6,2) * t59 - mrSges(6,3) * t62) * t56;
t74 = m(6) * (-t58 * t54 * qJ(5) + t35 * pkin(4) + qJDD(5) + ((2 * qJD(4)) + t44) * qJD(3) + t67) + t33 * t94 + qJD(3) * t45 + qJDD(3) * mrSges(6,2);
t69 = m(5) * (qJD(3) * t101 - t67) - t74;
t31 = (mrSges(5,2) * t62 - mrSges(5,3) * t59) * t56;
t87 = t31 + (-mrSges(4,1) * t62 + mrSges(4,2) * t59) * t56;
t89 = -mrSges(4,3) - mrSges(5,1);
t7 = m(4) * t76 + t87 * t94 + (-mrSges(4,2) + mrSges(5,3)) * qJDD(3) + (-t42 + t48) * qJD(3) + (mrSges(6,1) - t89) * t35 - t69;
t88 = -t62 * g(3) - t59 * t22;
t18 = -qJDD(3) * pkin(3) - t65 * qJ(4) + t30 * t95 + qJDD(4) - t88;
t81 = m(6) * (qJD(3) * t100 + (-t54 * t59 * t62 - qJDD(3)) * qJ(5) + (t34 - t78) * pkin(4) + t18) + t33 * t95 + t34 * mrSges(6,1);
t71 = m(5) * t18 + t81;
t83 = t46 - t47;
t90 = mrSges(5,2) - mrSges(6,3);
t8 = m(4) * t88 - t87 * t95 + t89 * t34 + (mrSges(4,1) - t90) * qJDD(3) + (t43 - t83) * qJD(3) - t71;
t98 = t59 * t7 + t62 * t8;
t96 = t34 * mrSges(6,2);
t93 = t62 * t47;
t4 = m(3) * t75 + t55 * mrSges(3,1) - t54 * mrSges(3,2) - t105;
t3 = m(3) * t86 - t54 * mrSges(3,1) - t55 * mrSges(3,2) - t59 * t8 + t62 * t7;
t2 = m(2) * t72 - t66 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t63 * t3 - t60 * t4;
t1 = m(2) * t77 + qJDD(1) * mrSges(2,1) - t66 * mrSges(2,2) + t60 * t3 + t63 * t4;
t5 = [-m(1) * g(1) - t61 * t1 + t64 * t2, t2, t3, t7, t35 * mrSges(5,2) - t96 + (t84 * t59 - t93) * t56 + t73, -t96 + (-t59 * t45 - t93) * t56 + t97; -m(1) * g(2) + t64 * t1 + t61 * t2, t1, t4, t8, -t31 * t94 - qJDD(3) * mrSges(5,3) - qJD(3) * t48 + (-mrSges(5,1) - mrSges(6,1)) * t35 + t69, -qJDD(3) * mrSges(6,3) - qJD(3) * t47 + t81; (-m(1) + t99) * g(3) + t98, t99 * g(3) + t98, -m(3) * g(3) + t98, t105, t34 * mrSges(5,1) + t83 * qJD(3) + t90 * qJDD(3) + t31 * t95 + t71, t35 * mrSges(6,1) + t74;];
f_new = t5;
