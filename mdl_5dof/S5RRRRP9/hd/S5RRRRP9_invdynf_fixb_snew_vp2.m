% Calculate vector of cutting forces with Newton-Euler
% S5RRRRP9
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRRP9_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP9_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP9_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP9_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP9_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP9_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP9_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP9_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP9_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:03:51
% EndTime: 2019-12-31 22:03:54
% DurationCPUTime: 1.20s
% Computational Cost: add. (13337->159), mult. (26408->202), div. (0->0), fcn. (17383->8), ass. (0->79)
t105 = cos(qJ(4));
t76 = sin(qJ(3));
t79 = cos(qJ(3));
t77 = sin(qJ(2));
t99 = qJD(1) * t77;
t61 = t79 * qJD(2) - t76 * t99;
t62 = t76 * qJD(2) + t79 * t99;
t75 = sin(qJ(4));
t45 = -t105 * t61 + t75 * t62;
t46 = t105 * t62 + t75 * t61;
t28 = t45 * mrSges(6,1) - t46 * mrSges(6,3);
t102 = -t45 * mrSges(5,1) - t46 * mrSges(5,2) - t28;
t104 = -mrSges(5,3) - mrSges(6,2);
t27 = t45 * pkin(4) - t46 * qJ(5);
t97 = qJD(1) * qJD(2);
t72 = t77 * t97;
t80 = cos(qJ(2));
t66 = t80 * qJDD(1) - t72;
t60 = qJDD(3) - t66;
t58 = qJDD(4) + t60;
t98 = t80 * qJD(1);
t71 = qJD(3) - t98;
t70 = qJD(4) + t71;
t69 = t70 ^ 2;
t93 = t80 * t97;
t65 = t77 * qJDD(1) + t93;
t44 = t61 * qJD(3) + t76 * qJDD(2) + t79 * t65;
t83 = qJD(1) ^ 2;
t78 = sin(qJ(1));
t81 = cos(qJ(1));
t95 = t78 * g(1) - t81 * g(2);
t56 = -qJDD(1) * pkin(1) - t83 * pkin(6) - t95;
t32 = (-t65 - t93) * pkin(7) + (-t66 + t72) * pkin(2) + t56;
t64 = (-pkin(2) * t80 - pkin(7) * t77) * qJD(1);
t82 = qJD(2) ^ 2;
t91 = -t81 * g(1) - t78 * g(2);
t57 = -t83 * pkin(1) + qJDD(1) * pkin(6) + t91;
t94 = -t77 * g(3) + t80 * t57;
t36 = -t82 * pkin(2) + qJDD(2) * pkin(7) + t64 * t98 + t94;
t92 = t79 * t32 - t76 * t36;
t16 = (t61 * t71 - t44) * pkin(8) + (t61 * t62 + t60) * pkin(3) + t92;
t101 = t76 * t32 + t79 * t36;
t43 = -t62 * qJD(3) + t79 * qJDD(2) - t76 * t65;
t50 = t71 * pkin(3) - t62 * pkin(8);
t59 = t61 ^ 2;
t18 = -t59 * pkin(3) + t43 * pkin(8) - t71 * t50 + t101;
t89 = t105 * t16 - t75 * t18;
t106 = m(6) * (-t58 * pkin(4) - t69 * qJ(5) + t46 * t27 + qJDD(5) - t89);
t24 = -t45 * qJD(4) + t105 * t44 + t75 * t43;
t37 = -t45 * mrSges(6,2) + t70 * mrSges(6,3);
t38 = -t70 * mrSges(5,2) - t45 * mrSges(5,3);
t10 = m(5) * t89 - t106 + (t38 + t37) * t70 + (mrSges(5,1) + mrSges(6,1)) * t58 + t102 * t46 + t104 * t24;
t47 = -t61 * mrSges(4,1) + t62 * mrSges(4,2);
t48 = -t71 * mrSges(4,2) + t61 * mrSges(4,3);
t103 = t105 * t18 + t75 * t16;
t23 = t46 * qJD(4) - t105 * t43 + t75 * t44;
t39 = t70 * mrSges(5,1) - t46 * mrSges(5,3);
t40 = -t70 * mrSges(6,1) + t46 * mrSges(6,2);
t96 = m(6) * (-t69 * pkin(4) + t58 * qJ(5) + 0.2e1 * qJD(5) * t70 - t45 * t27 + t103) + t70 * t40 + t58 * mrSges(6,3);
t9 = m(5) * t103 - t58 * mrSges(5,2) + t102 * t45 + t104 * t23 - t70 * t39 + t96;
t5 = m(4) * t92 + t60 * mrSges(4,1) - t44 * mrSges(4,3) + t105 * t10 - t62 * t47 + t71 * t48 + t75 * t9;
t49 = t71 * mrSges(4,1) - t62 * mrSges(4,3);
t6 = m(4) * t101 - t60 * mrSges(4,2) + t43 * mrSges(4,3) - t75 * t10 + t105 * t9 + t61 * t47 - t71 * t49;
t67 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t99;
t68 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t98;
t108 = m(3) * t56 - t66 * mrSges(3,1) + t65 * mrSges(3,2) + t79 * t5 + t76 * t6 + (t77 * t67 - t80 * t68) * qJD(1);
t63 = (-mrSges(3,1) * t80 + mrSges(3,2) * t77) * qJD(1);
t4 = m(3) * t94 - qJDD(2) * mrSges(3,2) + t66 * mrSges(3,3) - qJD(2) * t67 - t76 * t5 + t79 * t6 + t63 * t98;
t100 = -t80 * g(3) - t77 * t57;
t35 = -qJDD(2) * pkin(2) - t82 * pkin(7) + t64 * t99 - t100;
t86 = -t43 * pkin(3) - t59 * pkin(8) + t62 * t50 + t35;
t88 = t24 * mrSges(6,3) + t46 * t40 - m(6) * (-0.2e1 * qJD(5) * t46 + (t45 * t70 - t24) * qJ(5) + (t46 * t70 + t23) * pkin(4) + t86) - t23 * mrSges(6,1) - t45 * t37;
t85 = m(5) * t86 + t23 * mrSges(5,1) + t24 * mrSges(5,2) + t45 * t38 + t46 * t39 - t88;
t84 = m(4) * t35 - t43 * mrSges(4,1) + t44 * mrSges(4,2) - t61 * t48 + t62 * t49 + t85;
t8 = m(3) * t100 + qJDD(2) * mrSges(3,1) - t65 * mrSges(3,3) + qJD(2) * t68 - t63 * t99 - t84;
t107 = t77 * t4 + t80 * t8;
t2 = m(2) * t95 + qJDD(1) * mrSges(2,1) - t83 * mrSges(2,2) - t108;
t1 = m(2) * t91 - t83 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t80 * t4 - t77 * t8;
t3 = [-m(1) * g(1) + t81 * t1 - t78 * t2, t1, t4, t6, t9, -t23 * mrSges(6,2) - t45 * t28 + t96; -m(1) * g(2) + t78 * t1 + t81 * t2, t2, t8, t5, t10, -t88; (-m(1) - m(2)) * g(3) + t107, -m(2) * g(3) + t107, t108, t84, t85, -t58 * mrSges(6,1) + t24 * mrSges(6,2) + t46 * t28 - t70 * t37 + t106;];
f_new = t3;
