% Calculate vector of cutting forces with Newton-Euler
% S5RRRPR5
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRPR5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR5_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR5_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR5_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR5_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR5_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:13:06
% EndTime: 2019-12-31 21:13:11
% DurationCPUTime: 2.20s
% Computational Cost: add. (25866->167), mult. (57952->224), div. (0->0), fcn. (41067->10), ass. (0->86)
t103 = qJD(1) * qJD(2);
t82 = sin(qJ(2));
t86 = cos(qJ(2));
t68 = t82 * qJDD(1) + t86 * t103;
t69 = t86 * qJDD(1) - t82 * t103;
t105 = qJD(1) * t82;
t70 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t105;
t104 = qJD(1) * t86;
t71 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t104;
t88 = qJD(1) ^ 2;
t81 = sin(qJ(3));
t85 = cos(qJ(3));
t63 = (t81 * t86 + t82 * t85) * qJD(1);
t45 = -t63 * qJD(3) - t81 * t68 + t85 * t69;
t62 = (-t81 * t82 + t85 * t86) * qJD(1);
t46 = t62 * qJD(3) + t85 * t68 + t81 * t69;
t76 = qJD(2) + qJD(3);
t57 = -t76 * mrSges(4,2) + t62 * mrSges(4,3);
t59 = t76 * mrSges(4,1) - t63 * mrSges(4,3);
t72 = qJD(2) * pkin(2) - pkin(7) * t105;
t77 = t86 ^ 2;
t83 = sin(qJ(1));
t87 = cos(qJ(1));
t101 = t83 * g(1) - t87 * g(2);
t95 = -qJDD(1) * pkin(1) - t101;
t92 = -t69 * pkin(2) + t72 * t105 + (-pkin(7) * t77 - pkin(6)) * t88 + t95;
t110 = 2 * qJD(4);
t75 = qJDD(2) + qJDD(3);
t98 = -t87 * g(1) - t83 * g(2);
t65 = -t88 * pkin(1) + qJDD(1) * pkin(6) + t98;
t107 = t82 * t65;
t108 = pkin(2) * t88;
t39 = qJDD(2) * pkin(2) - t68 * pkin(7) - t107 + (pkin(7) * t103 + t82 * t108 - g(3)) * t86;
t100 = -t82 * g(3) + t86 * t65;
t40 = t69 * pkin(7) - qJD(2) * t72 - t77 * t108 + t100;
t99 = t85 * t39 - t81 * t40;
t19 = (t62 * t76 - t46) * qJ(4) + (t62 * t63 + t75) * pkin(3) + t99;
t106 = t81 * t39 + t85 * t40;
t58 = t76 * pkin(3) - t63 * qJ(4);
t61 = t62 ^ 2;
t21 = -t61 * pkin(3) + t45 * qJ(4) - t76 * t58 + t106;
t78 = sin(pkin(9));
t79 = cos(pkin(9));
t54 = t79 * t62 - t78 * t63;
t102 = t54 * t110 + t78 * t19 + t79 * t21;
t55 = t78 * t62 + t79 * t63;
t35 = -t54 * pkin(4) - t55 * pkin(8);
t74 = t76 ^ 2;
t16 = -t74 * pkin(4) + t75 * pkin(8) + t54 * t35 + t102;
t28 = t79 * t45 - t78 * t46;
t29 = t78 * t45 + t79 * t46;
t90 = -t45 * pkin(3) - t61 * qJ(4) + t63 * t58 + qJDD(4) + t92;
t17 = (-t54 * t76 - t29) * pkin(8) + (t55 * t76 - t28) * pkin(4) + t90;
t80 = sin(qJ(5));
t84 = cos(qJ(5));
t42 = -t80 * t55 + t84 * t76;
t23 = t42 * qJD(5) + t84 * t29 + t80 * t75;
t27 = qJDD(5) - t28;
t43 = t84 * t55 + t80 * t76;
t30 = -t42 * mrSges(6,1) + t43 * mrSges(6,2);
t51 = qJD(5) - t54;
t31 = -t51 * mrSges(6,2) + t42 * mrSges(6,3);
t13 = m(6) * (-t80 * t16 + t84 * t17) - t23 * mrSges(6,3) + t27 * mrSges(6,1) - t43 * t30 + t51 * t31;
t22 = -t43 * qJD(5) - t80 * t29 + t84 * t75;
t32 = t51 * mrSges(6,1) - t43 * mrSges(6,3);
t14 = m(6) * (t84 * t16 + t80 * t17) + t22 * mrSges(6,3) - t27 * mrSges(6,2) + t42 * t30 - t51 * t32;
t48 = -t76 * mrSges(5,2) + t54 * mrSges(5,3);
t49 = t76 * mrSges(5,1) - t55 * mrSges(5,3);
t94 = -m(5) * t90 + t28 * mrSges(5,1) - t29 * mrSges(5,2) - t84 * t13 - t80 * t14 + t54 * t48 - t55 * t49;
t91 = -m(4) * t92 + t45 * mrSges(4,1) - t46 * mrSges(4,2) + t62 * t57 - t63 * t59 + t94;
t111 = (t82 * t70 - t86 * t71) * qJD(1) + m(3) * (-t88 * pkin(6) + t95) - t69 * mrSges(3,1) + t68 * mrSges(3,2) - t91;
t34 = -t54 * mrSges(5,1) + t55 * mrSges(5,2);
t97 = -t79 * t19 + t78 * t21;
t93 = m(6) * (-t75 * pkin(4) - t74 * pkin(8) + (t110 + t35) * t55 + t97) - t22 * mrSges(6,1) + t23 * mrSges(6,2) - t42 * t31 + t43 * t32;
t10 = m(5) * (-0.2e1 * qJD(4) * t55 - t97) - t29 * mrSges(5,3) + t75 * mrSges(5,1) - t55 * t34 + t76 * t48 - t93;
t56 = -t62 * mrSges(4,1) + t63 * mrSges(4,2);
t9 = m(5) * t102 - t75 * mrSges(5,2) + t28 * mrSges(5,3) - t80 * t13 + t84 * t14 + t54 * t34 - t76 * t49;
t6 = m(4) * t99 + t75 * mrSges(4,1) - t46 * mrSges(4,3) + t79 * t10 - t63 * t56 + t76 * t57 + t78 * t9;
t67 = (-mrSges(3,1) * t86 + mrSges(3,2) * t82) * qJD(1);
t7 = m(4) * t106 - t75 * mrSges(4,2) + t45 * mrSges(4,3) - t78 * t10 + t62 * t56 - t76 * t59 + t79 * t9;
t4 = m(3) * (-t86 * g(3) - t107) - t68 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t67 * t105 + qJD(2) * t71 + t81 * t7 + t85 * t6;
t5 = m(3) * t100 - qJDD(2) * mrSges(3,2) + t69 * mrSges(3,3) - qJD(2) * t70 + t67 * t104 - t81 * t6 + t85 * t7;
t109 = t86 * t4 + t82 * t5;
t8 = m(2) * t101 + qJDD(1) * mrSges(2,1) - t88 * mrSges(2,2) - t111;
t1 = m(2) * t98 - t88 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t82 * t4 + t86 * t5;
t2 = [-m(1) * g(1) + t87 * t1 - t83 * t8, t1, t5, t7, t9, t14; -m(1) * g(2) + t83 * t1 + t87 * t8, t8, t4, t6, t10, t13; (-m(1) - m(2)) * g(3) + t109, -m(2) * g(3) + t109, t111, -t91, -t94, t93;];
f_new = t2;
