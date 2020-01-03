% Calculate vector of cutting forces with Newton-Euler
% S5RPRPR16
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
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
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRPR16_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR16_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR16_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR16_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR16_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR16_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR16_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:38:45
% EndTime: 2019-12-31 18:38:47
% DurationCPUTime: 0.60s
% Computational Cost: add. (2898->140), mult. (5658->167), div. (0->0), fcn. (2644->6), ass. (0->71)
t64 = qJD(1) ^ 2;
t58 = sin(qJ(3));
t61 = cos(qJ(3));
t87 = qJD(1) * qJD(3);
t81 = t61 * t87;
t41 = t58 * qJDD(1) + t81;
t52 = t58 * t87;
t42 = t61 * qJDD(1) - t52;
t89 = qJD(1) * t58;
t43 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t89;
t88 = t61 * qJD(1);
t44 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t88;
t45 = mrSges(5,1) * t89 - qJD(3) * mrSges(5,3);
t46 = mrSges(5,1) * t88 + qJD(3) * mrSges(5,2);
t57 = sin(qJ(5));
t60 = cos(qJ(5));
t101 = -2 * qJD(4);
t59 = sin(qJ(1));
t62 = cos(qJ(1));
t79 = -t62 * g(1) - t59 * g(2);
t75 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t79;
t68 = pkin(3) * t81 + t88 * t101 + t75 + (-t42 + t52) * qJ(4);
t38 = (pkin(3) * t58 - qJ(4) * t61) * qJD(1);
t63 = qJD(3) ^ 2;
t83 = t59 * g(1) - t62 * g(2);
t72 = -t64 * qJ(2) + qJDD(2) - t83;
t99 = -pkin(1) - pkin(6);
t27 = t99 * qJDD(1) + t72;
t94 = t61 * t27;
t71 = -t63 * qJ(4) + t38 * t88 + qJDD(4) - t94;
t96 = pkin(7) * t64;
t98 = pkin(3) + pkin(7);
t12 = t42 * pkin(4) - t98 * qJDD(3) + (pkin(4) * t87 + t61 * t96 - g(3)) * t58 + t71;
t36 = -t57 * qJD(3) + t60 * t89;
t20 = t36 * qJD(5) + t60 * qJDD(3) + t57 * t41;
t37 = t60 * qJD(3) + t57 * t89;
t21 = -t36 * mrSges(6,1) + t37 * mrSges(6,2);
t50 = qJD(5) + t88;
t22 = -t50 * mrSges(6,2) + t36 * mrSges(6,3);
t35 = qJDD(5) + t42;
t47 = pkin(4) * t88 - qJD(3) * pkin(7);
t56 = t58 ^ 2;
t9 = -t47 * t88 + t98 * t41 + (-pkin(4) * t56 + t99) * t64 + t68;
t7 = m(6) * (t60 * t12 - t57 * t9) - t20 * mrSges(6,3) + t35 * mrSges(6,1) - t37 * t21 + t50 * t22;
t19 = -t37 * qJD(5) - t57 * qJDD(3) + t60 * t41;
t23 = t50 * mrSges(6,1) - t37 * mrSges(6,3);
t8 = m(6) * (t57 * t12 + t60 * t9) + t19 * mrSges(6,3) - t35 * mrSges(6,2) + t36 * t21 - t50 * t23;
t85 = t99 * t64;
t67 = -(t58 * t45 + t61 * t46) * qJD(1) - t57 * t7 + t60 * t8 + m(5) * (t41 * pkin(3) + t68 + t85) - t42 * mrSges(5,3);
t92 = mrSges(4,1) - mrSges(5,2);
t65 = t92 * t41 + m(4) * (t85 + t75) + t43 * t89 + t44 * t88 + t42 * mrSges(4,2) + t67;
t103 = m(3) * (t64 * pkin(1) - t75) - t65;
t100 = -m(2) - m(3);
t95 = t58 * g(3);
t93 = mrSges(2,1) - mrSges(3,2);
t91 = -mrSges(2,2) + mrSges(3,3);
t90 = -mrSges(4,3) - mrSges(5,1);
t73 = -m(5) * (-qJDD(3) * pkin(3) + t71 - t95) - t57 * t8 - t60 * t7;
t39 = (-mrSges(5,2) * t58 - mrSges(5,3) * t61) * qJD(1);
t80 = qJD(1) * (-t39 - (mrSges(4,1) * t58 + mrSges(4,2) * t61) * qJD(1));
t3 = m(4) * (t94 + t95) + t90 * t42 + t92 * qJDD(3) + (t43 - t45) * qJD(3) + t61 * t80 + t73;
t82 = -t61 * g(3) + t58 * t27;
t66 = -t63 * pkin(3) + qJDD(3) * qJ(4) - t38 * t89 + t82;
t70 = -t19 * mrSges(6,1) - t36 * t22 + m(6) * (-t56 * t96 - t41 * pkin(4) + ((2 * qJD(4)) + t47) * qJD(3) + t66) + t20 * mrSges(6,2) + t37 * t23;
t69 = -m(5) * (qJD(3) * t101 - t66) + t70;
t5 = m(4) * t82 + t90 * t41 + (-mrSges(4,2) + mrSges(5,3)) * qJDD(3) + (-t44 + t46) * qJD(3) + t58 * t80 + t69;
t84 = -t58 * t3 + t61 * t5;
t74 = -m(3) * (-qJDD(1) * pkin(1) + t72) - t61 * t3 - t58 * t5;
t2 = m(2) * t79 + t91 * qJDD(1) - t93 * t64 - t103;
t1 = m(2) * t83 + t93 * qJDD(1) + t91 * t64 + t74;
t4 = [-m(1) * g(1) - t59 * t1 + t62 * t2, t2, -m(3) * g(3) + t84, t5, -t41 * mrSges(5,2) + t67, t8; -m(1) * g(2) + t62 * t1 + t59 * t2, t1, -t64 * mrSges(3,2) - qJDD(1) * mrSges(3,3) + t103, t3, t41 * mrSges(5,1) - qJDD(3) * mrSges(5,3) - qJD(3) * t46 + t39 * t89 - t69, t7; (-m(1) + t100) * g(3) + t84, t100 * g(3) + t84, qJDD(1) * mrSges(3,2) - t64 * mrSges(3,3) - t74, t65, t42 * mrSges(5,1) + qJDD(3) * mrSges(5,2) + qJD(3) * t45 + t39 * t88 - t73, t70;];
f_new = t4;
