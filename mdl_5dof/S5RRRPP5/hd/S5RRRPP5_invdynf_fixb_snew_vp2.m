% Calculate vector of cutting forces with Newton-Euler
% S5RRRPP5
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
% Datum: 2019-12-31 20:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRPP5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP5_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP5_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP5_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP5_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP5_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:57:36
% EndTime: 2019-12-31 20:57:39
% DurationCPUTime: 0.78s
% Computational Cost: add. (5385->162), mult. (11212->192), div. (0->0), fcn. (6750->6), ass. (0->72)
t80 = cos(qJ(2));
t102 = qJD(1) * t80;
t78 = sin(qJ(2));
t103 = qJD(1) * t78;
t110 = cos(qJ(3));
t77 = sin(qJ(3));
t54 = -t110 * t102 + t77 * t103;
t75 = qJD(2) + qJD(3);
t109 = t54 * t75;
t101 = qJD(1) * qJD(2);
t61 = t78 * qJDD(1) + t80 * t101;
t62 = t80 * qJDD(1) - t78 * t101;
t33 = -t54 * qJD(3) + t110 * t61 + t77 * t62;
t65 = qJD(2) * pkin(2) - pkin(7) * t103;
t76 = t80 ^ 2;
t82 = qJD(1) ^ 2;
t79 = sin(qJ(1));
t81 = cos(qJ(1));
t104 = t79 * g(1) - t81 * g(2);
t93 = qJDD(1) * pkin(1) + t104;
t85 = -t62 * pkin(2) + t65 * t103 - (pkin(7) * t76 + pkin(6)) * t82 - t93;
t118 = (-t33 + t109) * qJ(4) + t85;
t55 = (t110 * t78 + t77 * t80) * qJD(1);
t50 = -t75 * mrSges(5,1) + t55 * mrSges(5,2);
t74 = qJDD(2) + qJDD(3);
t94 = -t81 * g(1) - t79 * g(2);
t57 = -t82 * pkin(1) + qJDD(1) * pkin(6) + t94;
t108 = t78 * t57;
t111 = pkin(2) * t82;
t21 = qJDD(2) * pkin(2) - t61 * pkin(7) - t108 + (pkin(7) * t101 + t78 * t111 - g(3)) * t80;
t97 = -t78 * g(3) + t80 * t57;
t22 = t62 * pkin(7) - qJD(2) * t65 - t76 * t111 + t97;
t106 = t110 * t22 + t77 * t21;
t113 = 2 * qJD(4);
t39 = t54 * pkin(3) - t55 * qJ(4);
t73 = t75 ^ 2;
t89 = -t73 * pkin(3) + t74 * qJ(4) + t75 * t113 - t54 * t39 + t106;
t116 = m(5) * t89 + t74 * mrSges(5,3) + t75 * t50;
t63 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t103;
t64 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t102;
t32 = t55 * qJD(3) - t110 * t62 + t77 * t61;
t46 = -t75 * mrSges(4,2) - t54 * mrSges(4,3);
t49 = t75 * mrSges(4,1) - t55 * mrSges(4,3);
t114 = -0.2e1 * t55;
t51 = -t54 * mrSges(5,2) + t75 * mrSges(5,3);
t45 = t75 * mrSges(6,2) + t54 * mrSges(6,3);
t47 = -t75 * pkin(4) - t55 * qJ(5);
t48 = -t75 * mrSges(6,1) - t55 * mrSges(6,3);
t53 = t54 ^ 2;
t95 = m(6) * (-t53 * qJ(5) + qJDD(5) + (-pkin(3) - pkin(4)) * t32 + (-pkin(3) * t75 + t113 + t47) * t55 - t118) + t33 * mrSges(6,2) - t32 * mrSges(6,1) + t55 * t48 - t54 * t45;
t86 = t33 * mrSges(5,3) + t55 * t50 - m(5) * (qJD(4) * t114 + (t55 * t75 + t32) * pkin(3) + t118) - t32 * mrSges(5,1) - t54 * t51 + t95;
t84 = m(4) * t85 + t32 * mrSges(4,1) + t33 * mrSges(4,2) + t54 * t46 + t55 * t49 - t86;
t115 = (t78 * t63 - t80 * t64) * qJD(1) + m(3) * (-t82 * pkin(6) - t93) - t62 * mrSges(3,1) + t61 * mrSges(3,2) + t84;
t60 = (-mrSges(3,1) * t80 + mrSges(3,2) * t78) * qJD(1);
t40 = t54 * mrSges(5,1) - t55 * mrSges(5,3);
t105 = -t54 * mrSges(4,1) - t55 * mrSges(4,2) - t40;
t107 = -mrSges(4,3) - mrSges(5,2);
t41 = -t54 * mrSges(6,1) + t55 * mrSges(6,2);
t92 = t110 * t21 - t77 * t22;
t17 = -t74 * pkin(3) - t73 * qJ(4) + t55 * t39 + qJDD(4) - t92;
t99 = t75 * t45 + t74 * mrSges(6,1) - m(6) * (qJD(5) * t114 + (-t33 - t109) * qJ(5) + (t54 * t55 - t74) * pkin(4) + t17);
t90 = m(5) * t17 - t99;
t7 = m(4) * t92 + (t46 + t51) * t75 + (mrSges(4,1) + mrSges(5,1)) * t74 + (t41 + t105) * t55 + (mrSges(6,3) + t107) * t33 - t90;
t100 = m(6) * (-t53 * pkin(4) + t32 * qJ(5) + 0.2e1 * qJD(5) * t54 + t75 * t47 + t89) + t54 * t41 + t32 * mrSges(6,3);
t8 = m(4) * t106 + (-t49 + t48) * t75 + (-mrSges(4,2) + mrSges(6,2)) * t74 + t105 * t54 + t107 * t32 + t100 + t116;
t4 = m(3) * (-t80 * g(3) - t108) - t61 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t60 * t103 + qJD(2) * t64 + t77 * t8 + t110 * t7;
t5 = m(3) * t97 - qJDD(2) * mrSges(3,2) + t62 * mrSges(3,3) - qJD(2) * t63 + t60 * t102 + t110 * t8 - t77 * t7;
t112 = t80 * t4 + t78 * t5;
t88 = t74 * mrSges(6,2) + t75 * t48 + t100;
t6 = m(2) * t104 + qJDD(1) * mrSges(2,1) - t82 * mrSges(2,2) - t115;
t1 = m(2) * t94 - t82 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t78 * t4 + t80 * t5;
t2 = [-m(1) * g(1) + t81 * t1 - t79 * t6, t1, t5, t8, -t32 * mrSges(5,2) - t54 * t40 + t116 + t88, t88; -m(1) * g(2) + t79 * t1 + t81 * t6, t6, t4, t7, -t86, -t33 * mrSges(6,3) - t55 * t41 - t99; (-m(1) - m(2)) * g(3) + t112, -m(2) * g(3) + t112, t115, t84, -t74 * mrSges(5,1) - t75 * t51 + (t40 - t41) * t55 + (mrSges(5,2) - mrSges(6,3)) * t33 + t90, t95;];
f_new = t2;
