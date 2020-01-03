% Calculate vector of cutting forces with Newton-Euler
% S5RRRPP4
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
% Datum: 2019-12-31 20:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRPP4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP4_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP4_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP4_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP4_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP4_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:54:55
% EndTime: 2019-12-31 20:54:58
% DurationCPUTime: 1.36s
% Computational Cost: add. (14711->162), mult. (32798->210), div. (0->0), fcn. (22034->8), ass. (0->77)
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
t76 = sin(qJ(3));
t79 = cos(qJ(3));
t59 = (t76 * t80 + t77 * t79) * qJD(1);
t39 = -t59 * qJD(3) - t76 * t64 + t79 * t65;
t58 = (-t76 * t77 + t79 * t80) * qJD(1);
t40 = t58 * qJD(3) + t79 * t64 + t76 * t65;
t73 = qJD(2) + qJD(3);
t52 = -t73 * mrSges(4,2) + t58 * mrSges(4,3);
t54 = t73 * mrSges(4,1) - t59 * mrSges(4,3);
t101 = cos(pkin(8));
t75 = sin(pkin(8));
t23 = -t101 * t39 + t75 * t40;
t24 = t101 * t40 + t75 * t39;
t49 = -t101 * t58 + t75 * t59;
t42 = -t73 * mrSges(5,2) - t49 * mrSges(5,3);
t50 = t101 * t59 + t75 * t58;
t43 = t73 * mrSges(5,1) - t50 * mrSges(5,3);
t53 = t73 * pkin(3) - t59 * qJ(4);
t57 = t58 ^ 2;
t68 = qJD(2) * pkin(2) - pkin(7) * t100;
t74 = t80 ^ 2;
t78 = sin(qJ(1));
t81 = cos(qJ(1));
t95 = t78 * g(1) - t81 * g(2);
t90 = -qJDD(1) * pkin(1) - t95;
t87 = -t65 * pkin(2) + t68 * t100 + (-pkin(7) * t74 - pkin(6)) * t82 + t90;
t85 = -t39 * pkin(3) - t57 * qJ(4) + t59 * t53 + qJDD(4) + t87;
t44 = -t73 * mrSges(6,1) + t50 * mrSges(6,2);
t45 = -t49 * mrSges(6,2) + t73 * mrSges(6,3);
t88 = t24 * mrSges(6,3) + t50 * t44 - m(6) * (-0.2e1 * qJD(5) * t50 + (t49 * t73 - t24) * qJ(5) + (t50 * t73 + t23) * pkin(4) + t85) - t23 * mrSges(6,1) - t49 * t45;
t86 = m(5) * t85 + t23 * mrSges(5,1) + t24 * mrSges(5,2) + t49 * t42 + t50 * t43 - t88;
t84 = -m(4) * t87 + t39 * mrSges(4,1) - t40 * mrSges(4,2) + t58 * t52 - t59 * t54 - t86;
t110 = (t77 * t66 - t80 * t67) * qJD(1) + m(3) * (-t82 * pkin(6) + t90) - t65 * mrSges(3,1) + t64 * mrSges(3,2) - t84;
t109 = -2 * qJD(4);
t92 = -t81 * g(1) - t78 * g(2);
t61 = -t82 * pkin(1) + qJDD(1) * pkin(6) + t92;
t105 = t77 * t61;
t28 = t49 * mrSges(6,1) - t50 * mrSges(6,3);
t103 = -t49 * mrSges(5,1) - t50 * mrSges(5,2) - t28;
t104 = -mrSges(5,3) - mrSges(6,2);
t27 = t49 * pkin(4) - t50 * qJ(5);
t71 = t73 ^ 2;
t72 = qJDD(2) + qJDD(3);
t106 = pkin(2) * t82;
t33 = qJDD(2) * pkin(2) - t64 * pkin(7) - t105 + (pkin(7) * t98 + t77 * t106 - g(3)) * t80;
t94 = -t77 * g(3) + t80 * t61;
t34 = t65 * pkin(7) - qJD(2) * t68 - t74 * t106 + t94;
t93 = t79 * t33 - t76 * t34;
t16 = (t58 * t73 - t40) * qJ(4) + (t58 * t59 + t72) * pkin(3) + t93;
t102 = t76 * t33 + t79 * t34;
t18 = -t57 * pkin(3) + t39 * qJ(4) - t73 * t53 + t102;
t89 = t101 * t16 - t75 * t18;
t107 = m(6) * (-t72 * pkin(4) - t71 * qJ(5) + qJDD(5) + ((2 * qJD(4)) + t27) * t50 - t89);
t10 = m(5) * t89 - t107 + (t42 + t45) * t73 + (mrSges(5,1) + mrSges(6,1)) * t72 + (m(5) * t109 + t103) * t50 + t104 * t24;
t51 = -t58 * mrSges(4,1) + t59 * mrSges(4,2);
t96 = t101 * t18 + t49 * t109 + t75 * t16;
t97 = m(6) * (-t71 * pkin(4) + t72 * qJ(5) + 0.2e1 * qJD(5) * t73 - t49 * t27 + t96) + t73 * t44 + t72 * mrSges(6,3);
t9 = m(5) * t96 - t72 * mrSges(5,2) + t103 * t49 + t104 * t23 - t73 * t43 + t97;
t6 = m(4) * t93 + t72 * mrSges(4,1) - t40 * mrSges(4,3) + t101 * t10 - t59 * t51 + t73 * t52 + t75 * t9;
t63 = (-mrSges(3,1) * t80 + mrSges(3,2) * t77) * qJD(1);
t7 = m(4) * t102 - t72 * mrSges(4,2) + t39 * mrSges(4,3) - t75 * t10 + t101 * t9 + t58 * t51 - t73 * t54;
t4 = m(3) * (-t80 * g(3) - t105) - t64 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t63 * t100 + qJD(2) * t67 + t76 * t7 + t79 * t6;
t5 = m(3) * t94 - qJDD(2) * mrSges(3,2) + t65 * mrSges(3,3) - qJD(2) * t66 - t76 * t6 + t63 * t99 + t79 * t7;
t108 = t80 * t4 + t77 * t5;
t8 = m(2) * t95 + qJDD(1) * mrSges(2,1) - t82 * mrSges(2,2) - t110;
t1 = m(2) * t92 - t82 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t77 * t4 + t80 * t5;
t2 = [-m(1) * g(1) + t81 * t1 - t78 * t8, t1, t5, t7, t9, -t23 * mrSges(6,2) - t49 * t28 + t97; -m(1) * g(2) + t78 * t1 + t81 * t8, t8, t4, t6, t10, -t88; (-m(1) - m(2)) * g(3) + t108, -m(2) * g(3) + t108, t110, -t84, t86, -t72 * mrSges(6,1) + t24 * mrSges(6,2) + t50 * t28 - t73 * t45 + t107;];
f_new = t2;
