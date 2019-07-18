% Calculate vector of cutting forces with Newton-Euler
% S5RRPRR1
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
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
% Datum: 2019-07-18 17:22
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRPRR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(4,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR1_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR1_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_invdynf_fixb_snew_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR1_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR1_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR1_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:20:47
% EndTime: 2019-07-18 17:20:50
% DurationCPUTime: 0.89s
% Computational Cost: add. (7417->158), mult. (15999->198), div. (0->0), fcn. (9692->8), ass. (0->77)
t80 = sin(qJ(2));
t102 = qJD(1) * qJD(2);
t84 = cos(qJ(2));
t99 = t84 * t102;
t57 = t80 * qJDD(1) + t99;
t58 = t84 * qJDD(1) - t80 * t102;
t104 = qJD(1) * t80;
t62 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t104;
t103 = qJD(1) * t84;
t63 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t103;
t64 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t103;
t86 = qJD(1) ^ 2;
t111 = t84 ^ 2 * t86;
t61 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t104;
t107 = -pkin(3) - qJ(3);
t109 = t80 * t86;
t101 = qJD(1) * qJD(3);
t81 = sin(qJ(1));
t85 = cos(qJ(1));
t67 = -t85 * g(1) - t81 * g(2);
t110 = t80 * t67;
t92 = qJ(3) * t99 - 0.2e1 * t80 * t101 - t110 + (t109 * t84 + qJDD(2)) * pkin(1);
t23 = qJDD(2) * pkin(2) + t107 * t57 + (pkin(2) * t109 + pkin(3) * t102 - g(3)) * t84 + t92;
t115 = -pkin(1) - pkin(2);
t60 = qJD(2) * pkin(1) - qJ(3) * t104;
t65 = qJD(2) * pkin(2) - pkin(3) * t104;
t100 = -t80 * g(3) + t84 * t67;
t95 = t58 * qJ(3) + 0.2e1 * t84 * t101 + t100;
t24 = t58 * pkin(3) + t115 * t111 + (-t60 - t65) * qJD(2) + t95;
t79 = sin(qJ(4));
t83 = cos(qJ(4));
t106 = t79 * t23 + t83 * t24;
t46 = (t79 * t80 - t83 * t84) * qJD(1);
t47 = (t79 * t84 + t80 * t83) * qJD(1);
t73 = qJDD(2) + qJDD(4);
t16 = (t46 * t47 + t73) * pkin(4) + t106;
t33 = -t46 * qJD(4) + t83 * t57 + t79 * t58;
t74 = qJD(2) + qJD(4);
t66 = t81 * g(1) - t85 * g(2);
t94 = t60 * t104 + qJDD(3) - t66;
t87 = t65 * t104 + t107 * t111 + t115 * t58 + t94;
t18 = (t46 * t74 - t33) * pkin(4) + t87;
t78 = sin(qJ(5));
t82 = cos(qJ(5));
t40 = -t78 * t47 + t82 * t74;
t22 = t40 * qJD(5) + t82 * t33 + t78 * t73;
t41 = t82 * t47 + t78 * t74;
t28 = -t40 * mrSges(6,1) + t41 * mrSges(6,2);
t32 = -t47 * qJD(4) - t79 * t57 + t83 * t58;
t31 = qJDD(5) - t32;
t45 = qJD(5) + t46;
t34 = -t45 * mrSges(6,2) + t40 * mrSges(6,3);
t14 = m(6) * (-t78 * t16 + t82 * t18) - t22 * mrSges(6,3) + t31 * mrSges(6,1) - t41 * t28 + t45 * t34;
t21 = -t41 * qJD(5) - t78 * t33 + t82 * t73;
t35 = t45 * mrSges(6,1) - t41 * mrSges(6,3);
t15 = m(6) * (t82 * t16 + t78 * t18) + t21 * mrSges(6,3) - t31 * mrSges(6,2) + t40 * t28 - t45 * t35;
t42 = -t74 * mrSges(5,2) - t46 * mrSges(5,3);
t43 = t74 * mrSges(5,1) - t47 * mrSges(5,3);
t91 = -m(5) * t87 + t32 * mrSges(5,1) - t33 * mrSges(5,2) - t82 * t14 - t78 * t15 - t46 * t42 - t47 * t43;
t89 = -m(4) * (-t58 * pkin(1) - qJ(3) * t111 + t94) - t61 * t104 - t57 * mrSges(4,2) + t91;
t118 = (-t62 * t80 + (t63 + t64) * t84) * qJD(1) - t57 * mrSges(3,2) + (mrSges(3,1) + mrSges(4,1)) * t58 + t89;
t113 = t84 * g(3);
t55 = (-mrSges(4,1) * t84 + mrSges(4,2) * t80) * qJD(1);
t56 = (-mrSges(3,1) * t84 + mrSges(3,2) * t80) * qJD(1);
t39 = t46 * mrSges(5,1) + t47 * mrSges(5,2);
t97 = t83 * t23 - t79 * t24;
t90 = m(6) * ((-t47 ^ 2 - t74 ^ 2) * pkin(4) - t97) - t21 * mrSges(6,1) + t22 * mrSges(6,2) - t40 * t34 + t41 * t35;
t11 = m(5) * t97 + t73 * mrSges(5,1) - t33 * mrSges(5,3) - t47 * t39 + t74 * t42 - t90;
t9 = m(5) * t106 - t73 * mrSges(5,2) + t32 * mrSges(5,3) - t78 * t14 + t82 * t15 - t46 * t39 - t74 * t43;
t98 = t83 * t11 + t79 * t9 + m(4) * (-t57 * qJ(3) - t113 + t92) + qJD(2) * t63 + qJDD(2) * mrSges(4,1);
t4 = m(3) * (-t110 - t113) + qJDD(2) * mrSges(3,1) + qJD(2) * t64 + (-mrSges(3,3) - mrSges(4,3)) * t57 + (-t55 - t56) * t104 + t98;
t93 = -t79 * t11 + t83 * t9 + m(4) * (-pkin(1) * t111 - qJD(2) * t60 + t95) + t55 * t103 + t58 * mrSges(4,3);
t5 = m(3) * t100 + t58 * mrSges(3,3) + t56 * t103 + (-mrSges(3,2) - mrSges(4,2)) * qJDD(2) + (-t62 - t61) * qJD(2) + t93;
t114 = t84 * t4 + t80 * t5;
t6 = qJDD(1) * mrSges(2,1) + (m(2) + m(3)) * t66 - t86 * mrSges(2,2) + t118;
t1 = m(2) * t67 - t86 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t80 * t4 + t84 * t5;
t2 = [-m(1) * g(1) + t85 * t1 - t81 * t6, t1, t5, -qJDD(2) * mrSges(4,2) - qJD(2) * t61 + t93, t9, t15; -m(1) * g(2) + t81 * t1 + t85 * t6, t6, t4, -t57 * mrSges(4,3) - t55 * t104 + t98, t11, t14; (-m(1) - m(2)) * g(3) + t114, -m(2) * g(3) + t114, -m(3) * t66 - t118, -t58 * mrSges(4,1) - t63 * t103 - t89, -t91, t90;];
f_new  = t2;
