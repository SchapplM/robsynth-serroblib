% Calculate vector of cutting forces with Newton-Euler
% S5RRPPR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRPPR5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR5_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR5_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR5_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR5_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR5_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:28:56
% EndTime: 2019-12-31 19:28:59
% DurationCPUTime: 1.18s
% Computational Cost: add. (9657->164), mult. (22443->209), div. (0->0), fcn. (14264->8), ass. (0->79)
t85 = cos(qJ(2));
t108 = qJD(1) * t85;
t82 = sin(qJ(2));
t109 = qJD(1) * t82;
t110 = cos(pkin(8));
t80 = sin(pkin(8));
t59 = -t110 * t108 + t80 * t109;
t107 = qJD(2) * t59;
t106 = qJD(1) * qJD(2);
t66 = t82 * qJDD(1) + t85 * t106;
t67 = t85 * qJDD(1) - t82 * t106;
t48 = t110 * t66 + t80 * t67;
t68 = qJD(2) * pkin(2) - qJ(3) * t109;
t79 = t85 ^ 2;
t88 = qJD(1) ^ 2;
t83 = sin(qJ(1));
t86 = cos(qJ(1));
t111 = t83 * g(1) - t86 * g(2);
t99 = qJDD(1) * pkin(1) + t111;
t91 = -t67 * pkin(2) - (qJ(3) * t79 + pkin(6)) * t88 + t68 * t109 + qJDD(3) - t99;
t121 = (-t48 + t107) * qJ(4) + t91;
t69 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t109;
t70 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t108;
t47 = -t110 * t67 + t80 * t66;
t50 = -qJD(2) * mrSges(4,2) - t59 * mrSges(4,3);
t60 = (t110 * t82 + t80 * t85) * qJD(1);
t51 = qJD(2) * mrSges(4,1) - t60 * mrSges(4,3);
t117 = 2 * qJD(4);
t81 = sin(qJ(5));
t84 = cos(qJ(5));
t39 = t81 * t59 + t84 * t60;
t22 = -t39 * qJD(5) + t84 * t47 - t81 * t48;
t38 = t84 * t59 - t81 * t60;
t23 = t38 * qJD(5) + t81 * t47 + t84 * t48;
t76 = -qJD(2) + qJD(5);
t36 = -t76 * mrSges(6,2) + t38 * mrSges(6,3);
t37 = t76 * mrSges(6,1) - t39 * mrSges(6,3);
t54 = -qJD(2) * pkin(4) - t60 * pkin(7);
t58 = t59 ^ 2;
t101 = m(6) * (-t58 * pkin(7) + (-pkin(3) - pkin(4)) * t47 + (-pkin(3) * qJD(2) + t117 + t54) * t60 - t121) + t23 * mrSges(6,2) - t22 * mrSges(6,1) + t39 * t37 - t38 * t36;
t52 = -qJD(2) * mrSges(5,1) + t60 * mrSges(5,2);
t53 = -t59 * mrSges(5,2) + qJD(2) * mrSges(5,3);
t92 = t48 * mrSges(5,3) + t60 * t52 + t101 - m(5) * (-0.2e1 * qJD(4) * t60 + (qJD(2) * t60 + t47) * pkin(3) + t121) - t59 * t53 - t47 * mrSges(5,1);
t90 = m(4) * t91 + t47 * mrSges(4,1) + t48 * mrSges(4,2) + t59 * t50 + t60 * t51 - t92;
t119 = (t82 * t69 - t85 * t70) * qJD(1) + m(3) * (-t88 * pkin(6) - t99) - t67 * mrSges(3,1) + t66 * mrSges(3,2) + t90;
t118 = -2 * qJD(3);
t100 = -t86 * g(1) - t83 * g(2);
t63 = -t88 * pkin(1) + qJDD(1) * pkin(6) + t100;
t114 = t82 * t63;
t115 = pkin(2) * t88;
t30 = qJDD(2) * pkin(2) - t66 * qJ(3) - t114 + (qJ(3) * t106 + t82 * t115 - g(3)) * t85;
t104 = -t82 * g(3) + t85 * t63;
t31 = t67 * qJ(3) - qJD(2) * t68 - t79 * t115 + t104;
t105 = t110 * t31 + t59 * t118 + t80 * t30;
t43 = t59 * mrSges(5,1) - t60 * mrSges(5,3);
t112 = -t59 * mrSges(4,1) - t60 * mrSges(4,2) - t43;
t113 = -mrSges(4,3) - mrSges(5,2);
t42 = t59 * pkin(3) - t60 * qJ(4);
t87 = qJD(2) ^ 2;
t97 = t110 * t30 - t80 * t31;
t17 = -qJDD(2) * pkin(3) - t87 * qJ(4) + qJDD(4) - t97 + ((2 * qJD(3)) + t42) * t60;
t12 = (-t48 - t107) * pkin(7) + (t59 * t60 - qJDD(2)) * pkin(4) + t17;
t94 = -t87 * pkin(3) + qJDD(2) * qJ(4) + qJD(2) * t117 - t59 * t42 + t105;
t13 = -t58 * pkin(4) + t47 * pkin(7) + qJD(2) * t54 + t94;
t26 = -t38 * mrSges(6,1) + t39 * mrSges(6,2);
t75 = -qJDD(2) + qJDD(5);
t10 = m(6) * (t84 * t12 - t81 * t13) - t23 * mrSges(6,3) + t75 * mrSges(6,1) - t39 * t26 + t76 * t36;
t11 = m(6) * (t81 * t12 + t84 * t13) + t22 * mrSges(6,3) - t75 * mrSges(6,2) + t38 * t26 - t76 * t37;
t96 = m(5) * t94 + qJDD(2) * mrSges(5,3) + qJD(2) * t52 - t81 * t10 + t84 * t11;
t6 = m(4) * t105 - qJDD(2) * mrSges(4,2) - qJD(2) * t51 + t112 * t59 + t113 * t47 + t96;
t65 = (-mrSges(3,1) * t85 + mrSges(3,2) * t82) * qJD(1);
t95 = -m(5) * t17 - t84 * t10 - t81 * t11;
t7 = m(4) * t97 + (m(4) * t118 + t112) * t60 + t113 * t48 + (mrSges(4,1) + mrSges(5,1)) * qJDD(2) + (t50 + t53) * qJD(2) + t95;
t4 = m(3) * (-t85 * g(3) - t114) - t66 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t65 * t109 + qJD(2) * t70 + t80 * t6 + t110 * t7;
t5 = m(3) * t104 - qJDD(2) * mrSges(3,2) + t67 * mrSges(3,3) - qJD(2) * t69 + t65 * t108 + t110 * t6 - t80 * t7;
t116 = t85 * t4 + t82 * t5;
t8 = m(2) * t111 + qJDD(1) * mrSges(2,1) - t88 * mrSges(2,2) - t119;
t1 = m(2) * t100 - t88 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t82 * t4 + t85 * t5;
t2 = [-m(1) * g(1) + t86 * t1 - t83 * t8, t1, t5, t6, -t47 * mrSges(5,2) - t59 * t43 + t96, t11; -m(1) * g(2) + t83 * t1 + t86 * t8, t8, t4, t7, -t92, t10; (-m(1) - m(2)) * g(3) + t116, -m(2) * g(3) + t116, t119, t90, -qJDD(2) * mrSges(5,1) + t48 * mrSges(5,2) - qJD(2) * t53 + t60 * t43 - t95, t101;];
f_new = t2;
