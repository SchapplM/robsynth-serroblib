% Calculate vector of cutting forces with Newton-Euler
% S5RRRPR7
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
% Datum: 2019-12-31 21:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRPR7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR7_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR7_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR7_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR7_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR7_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR7_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR7_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:16:12
% EndTime: 2019-12-31 21:16:17
% DurationCPUTime: 2.11s
% Computational Cost: add. (25859->165), mult. (53924->223), div. (0->0), fcn. (37444->10), ass. (0->85)
t105 = qJD(1) * qJD(2);
t86 = sin(qJ(2));
t89 = cos(qJ(2));
t71 = t86 * qJDD(1) + t89 * t105;
t72 = t89 * qJDD(1) - t86 * t105;
t107 = qJD(1) * t86;
t73 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t107;
t106 = qJD(1) * t89;
t74 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t106;
t91 = qJD(1) ^ 2;
t110 = cos(qJ(3));
t85 = sin(qJ(3));
t65 = (t110 * t86 + t85 * t89) * qJD(1);
t46 = t65 * qJD(3) - t110 * t72 + t85 * t71;
t64 = -t110 * t106 + t85 * t107;
t47 = -t64 * qJD(3) + t110 * t71 + t85 * t72;
t80 = qJD(2) + qJD(3);
t75 = qJD(2) * pkin(2) - pkin(7) * t107;
t81 = t89 ^ 2;
t87 = sin(qJ(1));
t90 = cos(qJ(1));
t103 = t87 * g(1) - t90 * g(2);
t97 = -qJDD(1) * pkin(1) - t103;
t94 = -t72 * pkin(2) + t75 * t107 + (-pkin(7) * t81 - pkin(6)) * t91 + t97;
t24 = (t64 * t80 - t47) * qJ(4) + (t65 * t80 + t46) * pkin(3) + t94;
t100 = -t90 * g(1) - t87 * g(2);
t67 = -t91 * pkin(1) + qJDD(1) * pkin(6) + t100;
t109 = t86 * t67;
t111 = pkin(2) * t91;
t40 = qJDD(2) * pkin(2) - t71 * pkin(7) - t109 + (pkin(7) * t105 + t86 * t111 - g(3)) * t89;
t102 = -t86 * g(3) + t89 * t67;
t41 = t72 * pkin(7) - qJD(2) * t75 - t81 * t111 + t102;
t108 = t110 * t41 + t85 * t40;
t53 = t64 * pkin(3) - t65 * qJ(4);
t78 = t80 ^ 2;
t79 = qJDD(2) + qJDD(3);
t27 = -t78 * pkin(3) + t79 * qJ(4) - t64 * t53 + t108;
t82 = sin(pkin(9));
t83 = cos(pkin(9));
t59 = t83 * t65 + t82 * t80;
t101 = -0.2e1 * qJD(4) * t59 + t83 * t24 - t82 * t27;
t35 = t83 * t47 + t82 * t79;
t58 = -t82 * t65 + t83 * t80;
t15 = (t58 * t64 - t35) * pkin(8) + (t58 * t59 + t46) * pkin(4) + t101;
t104 = 0.2e1 * qJD(4) * t58 + t82 * t24 + t83 * t27;
t34 = -t82 * t47 + t83 * t79;
t51 = t64 * pkin(4) - t59 * pkin(8);
t57 = t58 ^ 2;
t16 = -t57 * pkin(4) + t34 * pkin(8) - t64 * t51 + t104;
t84 = sin(qJ(5));
t88 = cos(qJ(5));
t32 = t88 * t58 - t84 * t59;
t23 = t32 * qJD(5) + t84 * t34 + t88 * t35;
t33 = t84 * t58 + t88 * t59;
t29 = -t32 * mrSges(6,1) + t33 * mrSges(6,2);
t63 = qJD(5) + t64;
t30 = -t63 * mrSges(6,2) + t32 * mrSges(6,3);
t44 = qJDD(5) + t46;
t13 = m(6) * (t88 * t15 - t84 * t16) - t23 * mrSges(6,3) + t44 * mrSges(6,1) - t33 * t29 + t63 * t30;
t22 = -t33 * qJD(5) + t88 * t34 - t84 * t35;
t31 = t63 * mrSges(6,1) - t33 * mrSges(6,3);
t14 = m(6) * (t84 * t15 + t88 * t16) + t22 * mrSges(6,3) - t44 * mrSges(6,2) + t32 * t29 - t63 * t31;
t37 = -t58 * mrSges(5,1) + t59 * mrSges(5,2);
t49 = -t64 * mrSges(5,2) + t58 * mrSges(5,3);
t10 = m(5) * t101 + t46 * mrSges(5,1) - t35 * mrSges(5,3) + t88 * t13 + t84 * t14 - t59 * t37 + t64 * t49;
t50 = t64 * mrSges(5,1) - t59 * mrSges(5,3);
t11 = m(5) * t104 - t46 * mrSges(5,2) + t34 * mrSges(5,3) - t84 * t13 + t88 * t14 + t58 * t37 - t64 * t50;
t60 = -t80 * mrSges(4,2) - t64 * mrSges(4,3);
t61 = t80 * mrSges(4,1) - t65 * mrSges(4,3);
t95 = m(4) * t94 + t46 * mrSges(4,1) + t47 * mrSges(4,2) + t83 * t10 + t82 * t11 + t64 * t60 + t65 * t61;
t113 = (t86 * t73 - t89 * t74) * qJD(1) + m(3) * (-t91 * pkin(6) + t97) - t72 * mrSges(3,1) + t71 * mrSges(3,2) + t95;
t54 = t64 * mrSges(4,1) + t65 * mrSges(4,2);
t99 = t110 * t40 - t85 * t41;
t26 = -t79 * pkin(3) - t78 * qJ(4) + t65 * t53 + qJDD(4) - t99;
t96 = t22 * mrSges(6,1) + t32 * t30 - m(6) * (-t34 * pkin(4) - t57 * pkin(8) + t59 * t51 + t26) - t23 * mrSges(6,2) - t33 * t31;
t92 = m(5) * t26 - t34 * mrSges(5,1) + t35 * mrSges(5,2) - t58 * t49 + t59 * t50 - t96;
t12 = m(4) * t99 + t79 * mrSges(4,1) - t47 * mrSges(4,3) - t65 * t54 + t80 * t60 - t92;
t7 = m(4) * t108 - t79 * mrSges(4,2) - t46 * mrSges(4,3) - t82 * t10 + t83 * t11 - t64 * t54 - t80 * t61;
t70 = (-mrSges(3,1) * t89 + mrSges(3,2) * t86) * qJD(1);
t4 = m(3) * (-t89 * g(3) - t109) - t71 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t70 * t107 + qJD(2) * t74 + t85 * t7 + t110 * t12;
t5 = m(3) * t102 - qJDD(2) * mrSges(3,2) + t72 * mrSges(3,3) - qJD(2) * t73 + t70 * t106 + t110 * t7 - t85 * t12;
t112 = t89 * t4 + t86 * t5;
t6 = m(2) * t103 + qJDD(1) * mrSges(2,1) - t91 * mrSges(2,2) - t113;
t1 = m(2) * t100 - t91 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t86 * t4 + t89 * t5;
t2 = [-m(1) * g(1) + t90 * t1 - t87 * t6, t1, t5, t7, t11, t14; -m(1) * g(2) + t87 * t1 + t90 * t6, t6, t4, t12, t10, t13; (-m(1) - m(2)) * g(3) + t112, -m(2) * g(3) + t112, t113, t95, t92, -t96;];
f_new = t2;
