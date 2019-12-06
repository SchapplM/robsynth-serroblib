% Calculate vector of cutting forces with Newton-Euler
% S5RRRPR1
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
% Datum: 2019-12-05 18:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRPR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR1_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR1_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR1_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR1_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR1_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:37:55
% EndTime: 2019-12-05 18:38:01
% DurationCPUTime: 3.03s
% Computational Cost: add. (37154->167), mult. (85098->224), div. (0->0), fcn. (60296->10), ass. (0->85)
t106 = qJD(1) * qJD(2);
t85 = sin(qJ(2));
t89 = cos(qJ(2));
t70 = t85 * qJDD(1) + t89 * t106;
t71 = t89 * qJDD(1) - t85 * t106;
t108 = qJD(1) * t85;
t72 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t108;
t107 = qJD(1) * t89;
t73 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t107;
t91 = qJD(1) ^ 2;
t84 = sin(qJ(3));
t88 = cos(qJ(3));
t65 = (t84 * t89 + t85 * t88) * qJD(1);
t45 = -t65 * qJD(3) - t84 * t70 + t88 * t71;
t64 = (-t84 * t85 + t88 * t89) * qJD(1);
t46 = t64 * qJD(3) + t88 * t70 + t84 * t71;
t79 = qJD(2) + qJD(3);
t59 = -t79 * mrSges(4,2) + t64 * mrSges(4,3);
t61 = t79 * mrSges(4,1) - t65 * mrSges(4,3);
t81 = sin(pkin(9));
t82 = cos(pkin(9));
t30 = t82 * t45 - t81 * t46;
t31 = t81 * t45 + t82 * t46;
t56 = t82 * t64 - t81 * t65;
t48 = -t79 * mrSges(5,2) + t56 * mrSges(5,3);
t57 = t81 * t64 + t82 * t65;
t49 = t79 * mrSges(5,1) - t57 * mrSges(5,3);
t60 = t79 * pkin(3) - t65 * qJ(4);
t63 = t64 ^ 2;
t74 = qJD(2) * pkin(2) - pkin(7) * t108;
t80 = t89 ^ 2;
t86 = sin(qJ(1));
t90 = cos(qJ(1));
t104 = t86 * g(1) - t90 * g(2);
t98 = -qJDD(1) * pkin(1) - t104;
t96 = -t71 * pkin(2) + t74 * t108 + (-pkin(7) * t80 - pkin(6)) * t91 + t98;
t94 = -t45 * pkin(3) - t63 * qJ(4) + t65 * t60 + qJDD(4) + t96;
t83 = sin(qJ(5));
t87 = cos(qJ(5));
t36 = t83 * t56 + t87 * t57;
t18 = -t36 * qJD(5) + t87 * t30 - t83 * t31;
t35 = t87 * t56 - t83 * t57;
t19 = t35 * qJD(5) + t83 * t30 + t87 * t31;
t76 = qJD(5) + t79;
t32 = -t76 * mrSges(6,2) + t35 * mrSges(6,3);
t33 = t76 * mrSges(6,1) - t36 * mrSges(6,3);
t50 = t79 * pkin(4) - t57 * pkin(8);
t53 = t56 ^ 2;
t97 = t18 * mrSges(6,1) + t35 * t32 - m(6) * (-t30 * pkin(4) - t53 * pkin(8) + t57 * t50 + t94) - t19 * mrSges(6,2) - t36 * t33;
t95 = -m(5) * t94 + t30 * mrSges(5,1) - t31 * mrSges(5,2) + t56 * t48 - t57 * t49 + t97;
t93 = -m(4) * t96 + t45 * mrSges(4,1) - t46 * mrSges(4,2) + t64 * t59 - t65 * t61 + t95;
t113 = (t85 * t72 - t89 * t73) * qJD(1) + m(3) * (-t91 * pkin(6) + t98) - t71 * mrSges(3,1) + t70 * mrSges(3,2) - t93;
t100 = -t90 * g(1) - t86 * g(2);
t67 = -t91 * pkin(1) + qJDD(1) * pkin(6) + t100;
t110 = t85 * t67;
t111 = pkin(2) * t91;
t41 = qJDD(2) * pkin(2) - t70 * pkin(7) - t110 + (pkin(7) * t106 + t85 * t111 - g(3)) * t89;
t103 = -t85 * g(3) + t89 * t67;
t42 = t71 * pkin(7) - qJD(2) * t74 - t80 * t111 + t103;
t102 = t88 * t41 - t84 * t42;
t58 = -t64 * mrSges(4,1) + t65 * mrSges(4,2);
t78 = qJDD(2) + qJDD(3);
t22 = (t64 * t79 - t46) * qJ(4) + (t64 * t65 + t78) * pkin(3) + t102;
t109 = t84 * t41 + t88 * t42;
t24 = -t63 * pkin(3) + t45 * qJ(4) - t79 * t60 + t109;
t101 = -0.2e1 * qJD(4) * t57 + t82 * t22 - t81 * t24;
t13 = (t56 * t79 - t31) * pkin(8) + (t56 * t57 + t78) * pkin(4) + t101;
t105 = 0.2e1 * qJD(4) * t56 + t81 * t22 + t82 * t24;
t14 = -t53 * pkin(4) + t30 * pkin(8) - t79 * t50 + t105;
t26 = -t35 * mrSges(6,1) + t36 * mrSges(6,2);
t75 = qJDD(5) + t78;
t11 = m(6) * (t87 * t13 - t83 * t14) - t19 * mrSges(6,3) + t75 * mrSges(6,1) - t36 * t26 + t76 * t32;
t12 = m(6) * (t83 * t13 + t87 * t14) + t18 * mrSges(6,3) - t75 * mrSges(6,2) + t35 * t26 - t76 * t33;
t37 = -t56 * mrSges(5,1) + t57 * mrSges(5,2);
t8 = m(5) * t101 + t78 * mrSges(5,1) - t31 * mrSges(5,3) + t87 * t11 + t83 * t12 - t57 * t37 + t79 * t48;
t9 = m(5) * t105 - t78 * mrSges(5,2) + t30 * mrSges(5,3) - t83 * t11 + t87 * t12 + t56 * t37 - t79 * t49;
t6 = m(4) * t102 + t78 * mrSges(4,1) - t46 * mrSges(4,3) - t65 * t58 + t79 * t59 + t82 * t8 + t81 * t9;
t69 = (-mrSges(3,1) * t89 + mrSges(3,2) * t85) * qJD(1);
t7 = m(4) * t109 - t78 * mrSges(4,2) + t45 * mrSges(4,3) + t64 * t58 - t79 * t61 - t81 * t8 + t82 * t9;
t4 = m(3) * (-t89 * g(3) - t110) - t70 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t69 * t108 + qJD(2) * t73 + t84 * t7 + t88 * t6;
t5 = m(3) * t103 - qJDD(2) * mrSges(3,2) + t71 * mrSges(3,3) - qJD(2) * t72 + t69 * t107 - t84 * t6 + t88 * t7;
t112 = t89 * t4 + t85 * t5;
t10 = m(2) * t104 + qJDD(1) * mrSges(2,1) - t91 * mrSges(2,2) - t113;
t1 = m(2) * t100 - t91 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t85 * t4 + t89 * t5;
t2 = [-m(1) * g(1) + t90 * t1 - t86 * t10, t1, t5, t7, t9, t12; -m(1) * g(2) + t86 * t1 + t90 * t10, t10, t4, t6, t8, t11; (-m(1) - m(2)) * g(3) + t112, -m(2) * g(3) + t112, t113, -t93, -t95, -t97;];
f_new = t2;
