% Calculate vector of cutting forces with Newton-Euler
% S5RRRRP8
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
% Datum: 2019-12-31 22:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRRP8_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP8_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP8_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP8_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP8_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP8_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP8_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:59:43
% EndTime: 2019-12-31 21:59:46
% DurationCPUTime: 1.36s
% Computational Cost: add. (13826->160), mult. (27680->202), div. (0->0), fcn. (18454->8), ass. (0->78)
t78 = sin(qJ(2));
t104 = qJD(1) * t78;
t77 = sin(qJ(3));
t81 = cos(qJ(3));
t64 = t77 * qJD(2) + t81 * t104;
t102 = qJD(1) * qJD(2);
t82 = cos(qJ(2));
t94 = t82 * t102;
t67 = t78 * qJDD(1) + t94;
t45 = -t64 * qJD(3) + t81 * qJDD(2) - t77 * t67;
t63 = t81 * qJD(2) - t77 * t104;
t46 = t63 * qJD(3) + t77 * qJDD(2) + t81 * t67;
t76 = sin(qJ(4));
t80 = cos(qJ(4));
t49 = t76 * t63 + t80 * t64;
t25 = -t49 * qJD(4) + t80 * t45 - t76 * t46;
t48 = t80 * t63 - t76 * t64;
t26 = t48 * qJD(4) + t76 * t45 + t80 * t46;
t103 = t82 * qJD(1);
t72 = qJD(3) - t103;
t71 = qJD(4) + t72;
t39 = -t71 * mrSges(6,2) + t48 * mrSges(6,3);
t40 = -t71 * mrSges(5,2) + t48 * mrSges(5,3);
t43 = t71 * mrSges(5,1) - t49 * mrSges(5,3);
t85 = qJD(1) ^ 2;
t79 = sin(qJ(1));
t83 = cos(qJ(1));
t91 = -t83 * g(1) - t79 * g(2);
t59 = -t85 * pkin(1) + qJDD(1) * pkin(6) + t91;
t105 = -t82 * g(3) - t78 * t59;
t66 = (-pkin(2) * t82 - pkin(7) * t78) * qJD(1);
t84 = qJD(2) ^ 2;
t37 = -qJDD(2) * pkin(2) - t84 * pkin(7) + t66 * t104 - t105;
t53 = t72 * pkin(3) - t64 * pkin(8);
t61 = t63 ^ 2;
t87 = -t45 * pkin(3) - t61 * pkin(8) + t64 * t53 + t37;
t41 = t71 * pkin(4) - t49 * qJ(5);
t42 = t71 * mrSges(6,1) - t49 * mrSges(6,3);
t47 = t48 ^ 2;
t99 = m(6) * (-t25 * pkin(4) - t47 * qJ(5) + t49 * t41 + qJDD(5) + t87) + t26 * mrSges(6,2) + t49 * t42;
t114 = m(5) * t87 + t26 * mrSges(5,2) + t49 * t43 - (t40 + t39) * t48 - (mrSges(5,1) + mrSges(6,1)) * t25 + t99;
t51 = -t72 * mrSges(4,2) + t63 * mrSges(4,3);
t52 = t72 * mrSges(4,1) - t64 * mrSges(4,3);
t113 = m(4) * t37 - t45 * mrSges(4,1) + t46 * mrSges(4,2) - t63 * t51 + t64 * t52 + t114;
t50 = -t63 * mrSges(4,1) + t64 * mrSges(4,2);
t73 = t78 * t102;
t68 = t82 * qJDD(1) - t73;
t62 = qJDD(3) - t68;
t60 = qJDD(4) + t62;
t98 = t79 * g(1) - t83 * g(2);
t58 = -qJDD(1) * pkin(1) - t85 * pkin(6) - t98;
t34 = (-t67 - t94) * pkin(7) + (-t68 + t73) * pkin(2) + t58;
t97 = -t78 * g(3) + t82 * t59;
t38 = -t84 * pkin(2) + qJDD(2) * pkin(7) + t66 * t103 + t97;
t92 = t81 * t34 - t77 * t38;
t17 = (t63 * t72 - t46) * pkin(8) + (t63 * t64 + t62) * pkin(3) + t92;
t107 = t77 * t34 + t81 * t38;
t19 = -t61 * pkin(3) + t45 * pkin(8) - t72 * t53 + t107;
t93 = t80 * t17 - t76 * t19;
t101 = m(6) * (-0.2e1 * qJD(5) * t49 + (t48 * t71 - t26) * qJ(5) + (t48 * t49 + t60) * pkin(4) + t93) + t71 * t39 + t60 * mrSges(6,1);
t30 = -t48 * mrSges(6,1) + t49 * mrSges(6,2);
t31 = -t48 * mrSges(5,1) + t49 * mrSges(5,2);
t7 = m(5) * t93 + t60 * mrSges(5,1) + t71 * t40 + (-t31 - t30) * t49 + (-mrSges(5,3) - mrSges(6,3)) * t26 + t101;
t108 = t76 * t17 + t80 * t19;
t100 = m(6) * (-t47 * pkin(4) + t25 * qJ(5) + 0.2e1 * qJD(5) * t48 - t71 * t41 + t108) + t25 * mrSges(6,3) + t48 * t30;
t8 = m(5) * t108 + t25 * mrSges(5,3) + t48 * t31 + (-t43 - t42) * t71 + (-mrSges(5,2) - mrSges(6,2)) * t60 + t100;
t5 = m(4) * t92 + t62 * mrSges(4,1) - t46 * mrSges(4,3) - t64 * t50 + t72 * t51 + t80 * t7 + t76 * t8;
t6 = m(4) * t107 - t62 * mrSges(4,2) + t45 * mrSges(4,3) + t63 * t50 - t72 * t52 - t76 * t7 + t80 * t8;
t69 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t104;
t70 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t103;
t111 = m(3) * t58 - t68 * mrSges(3,1) + t67 * mrSges(3,2) + t81 * t5 + t77 * t6 + (t78 * t69 - t82 * t70) * qJD(1);
t65 = (-mrSges(3,1) * t82 + mrSges(3,2) * t78) * qJD(1);
t10 = m(3) * t105 + qJDD(2) * mrSges(3,1) - t67 * mrSges(3,3) + qJD(2) * t70 - t65 * t104 - t113;
t4 = m(3) * t97 - qJDD(2) * mrSges(3,2) + t68 * mrSges(3,3) - qJD(2) * t69 + t65 * t103 - t77 * t5 + t81 * t6;
t110 = t82 * t10 + t78 * t4;
t2 = m(2) * t98 + qJDD(1) * mrSges(2,1) - t85 * mrSges(2,2) - t111;
t1 = m(2) * t91 - t85 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t78 * t10 + t82 * t4;
t3 = [-m(1) * g(1) + t83 * t1 - t79 * t2, t1, t4, t6, t8, -t60 * mrSges(6,2) - t71 * t42 + t100; -m(1) * g(2) + t79 * t1 + t83 * t2, t2, t10, t5, t7, -t26 * mrSges(6,3) - t49 * t30 + t101; (-m(1) - m(2)) * g(3) + t110, -m(2) * g(3) + t110, t111, t113, t114, -t25 * mrSges(6,1) - t48 * t39 + t99;];
f_new = t3;
