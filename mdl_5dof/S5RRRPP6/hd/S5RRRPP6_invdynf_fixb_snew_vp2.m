% Calculate vector of cutting forces with Newton-Euler
% S5RRRPP6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Datum: 2019-12-31 21:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRPP6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP6_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP6_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP6_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP6_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP6_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:00:02
% EndTime: 2019-12-31 21:00:05
% DurationCPUTime: 1.19s
% Computational Cost: add. (12692->159), mult. (26034->204), div. (0->0), fcn. (16967->8), ass. (0->78)
t100 = cos(pkin(8));
t74 = sin(qJ(3));
t77 = cos(qJ(3));
t75 = sin(qJ(2));
t99 = qJD(1) * t75;
t61 = t77 * qJD(2) - t74 * t99;
t78 = cos(qJ(2));
t97 = qJD(1) * qJD(2);
t91 = t78 * t97;
t65 = t75 * qJDD(1) + t91;
t45 = t61 * qJD(3) + t74 * qJDD(2) + t77 * t65;
t62 = t74 * qJD(2) + t77 * t99;
t48 = -t61 * mrSges(4,1) + t62 * mrSges(4,2);
t98 = t78 * qJD(1);
t70 = qJD(3) - t98;
t49 = -t70 * mrSges(4,2) + t61 * mrSges(4,3);
t92 = t75 * t97;
t66 = t78 * qJDD(1) - t92;
t60 = qJDD(3) - t66;
t73 = sin(pkin(8));
t46 = -t100 * t61 + t73 * t62;
t47 = t100 * t62 + t73 * t61;
t28 = t46 * mrSges(6,1) - t47 * mrSges(6,3);
t103 = -t46 * mrSges(5,1) - t47 * mrSges(5,2) - t28;
t104 = -mrSges(5,3) - mrSges(6,2);
t44 = -t62 * qJD(3) + t77 * qJDD(2) - t74 * t65;
t23 = -t100 * t44 + t73 * t45;
t38 = t70 * mrSges(5,1) - t47 * mrSges(5,3);
t107 = -2 * qJD(4);
t81 = qJD(1) ^ 2;
t76 = sin(qJ(1));
t79 = cos(qJ(1));
t94 = t76 * g(1) - t79 * g(2);
t56 = -qJDD(1) * pkin(1) - t81 * pkin(6) - t94;
t32 = (-t65 - t91) * pkin(7) + (-t66 + t92) * pkin(2) + t56;
t64 = (-pkin(2) * t78 - pkin(7) * t75) * qJD(1);
t80 = qJD(2) ^ 2;
t89 = -t79 * g(1) - t76 * g(2);
t57 = -t81 * pkin(1) + qJDD(1) * pkin(6) + t89;
t93 = -t75 * g(3) + t78 * t57;
t36 = -t80 * pkin(2) + qJDD(2) * pkin(7) + t64 * t98 + t93;
t90 = t77 * t32 - t74 * t36;
t16 = (t61 * t70 - t45) * qJ(4) + (t61 * t62 + t60) * pkin(3) + t90;
t102 = t74 * t32 + t77 * t36;
t50 = t70 * pkin(3) - t62 * qJ(4);
t59 = t61 ^ 2;
t18 = -t59 * pkin(3) + t44 * qJ(4) - t70 * t50 + t102;
t95 = t100 * t18 + t46 * t107 + t73 * t16;
t27 = t46 * pkin(4) - t47 * qJ(5);
t39 = -t70 * mrSges(6,1) + t47 * mrSges(6,2);
t69 = t70 ^ 2;
t96 = m(6) * (-t69 * pkin(4) + t60 * qJ(5) + 0.2e1 * qJD(5) * t70 - t46 * t27 + t95) + t70 * t39 + t60 * mrSges(6,3);
t7 = m(5) * t95 - t60 * mrSges(5,2) + t103 * t46 + t104 * t23 - t70 * t38 + t96;
t87 = t100 * t16 - t73 * t18;
t105 = m(6) * (-t60 * pkin(4) - t69 * qJ(5) + qJDD(5) + ((2 * qJD(4)) + t27) * t47 - t87);
t24 = t100 * t45 + t73 * t44;
t37 = -t70 * mrSges(5,2) - t46 * mrSges(5,3);
t40 = -t46 * mrSges(6,2) + t70 * mrSges(6,3);
t8 = m(5) * t87 - t105 + (t37 + t40) * t70 + (mrSges(5,1) + mrSges(6,1)) * t60 + (m(5) * t107 + t103) * t47 + t104 * t24;
t5 = m(4) * t90 + t60 * mrSges(4,1) - t45 * mrSges(4,3) + t100 * t8 - t62 * t48 + t70 * t49 + t73 * t7;
t51 = t70 * mrSges(4,1) - t62 * mrSges(4,3);
t6 = m(4) * t102 - t60 * mrSges(4,2) + t44 * mrSges(4,3) + t100 * t7 + t61 * t48 - t70 * t51 - t73 * t8;
t67 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t99;
t68 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t98;
t108 = m(3) * t56 - t66 * mrSges(3,1) + t65 * mrSges(3,2) + t77 * t5 + t74 * t6 + (t75 * t67 - t78 * t68) * qJD(1);
t101 = -t78 * g(3) - t75 * t57;
t63 = (-mrSges(3,1) * t78 + mrSges(3,2) * t75) * qJD(1);
t35 = -qJDD(2) * pkin(2) - t80 * pkin(7) + t64 * t99 - t101;
t84 = -t44 * pkin(3) - t59 * qJ(4) + t62 * t50 + qJDD(4) + t35;
t86 = t24 * mrSges(6,3) + t47 * t39 - m(6) * (-0.2e1 * qJD(5) * t47 + (t46 * t70 - t24) * qJ(5) + (t47 * t70 + t23) * pkin(4) + t84) - t23 * mrSges(6,1) - t46 * t40;
t83 = m(5) * t84 + t23 * mrSges(5,1) + t24 * mrSges(5,2) + t46 * t37 + t47 * t38 - t86;
t82 = m(4) * t35 - t44 * mrSges(4,1) + t45 * mrSges(4,2) - t61 * t49 + t62 * t51 + t83;
t10 = m(3) * t101 + qJDD(2) * mrSges(3,1) - t65 * mrSges(3,3) + qJD(2) * t68 - t63 * t99 - t82;
t4 = m(3) * t93 - qJDD(2) * mrSges(3,2) + t66 * mrSges(3,3) - qJD(2) * t67 - t74 * t5 + t77 * t6 + t63 * t98;
t106 = t78 * t10 + t75 * t4;
t2 = m(2) * t94 + qJDD(1) * mrSges(2,1) - t81 * mrSges(2,2) - t108;
t1 = m(2) * t89 - t81 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t75 * t10 + t78 * t4;
t3 = [-m(1) * g(1) + t79 * t1 - t76 * t2, t1, t4, t6, t7, -t23 * mrSges(6,2) - t46 * t28 + t96; -m(1) * g(2) + t76 * t1 + t79 * t2, t2, t10, t5, t8, -t86; (-m(1) - m(2)) * g(3) + t106, -m(2) * g(3) + t106, t108, t82, t83, -t60 * mrSges(6,1) + t24 * mrSges(6,2) + t47 * t28 - t70 * t40 + t105;];
f_new = t3;
