% Calculate vector of cutting forces with Newton-Euler
% S5RRPRR13
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRPRR13_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR13_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR13_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR13_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR13_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR13_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR13_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR13_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR13_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:32:14
% EndTime: 2019-12-31 20:32:20
% DurationCPUTime: 2.37s
% Computational Cost: add. (28963->164), mult. (63052->218), div. (0->0), fcn. (43477->10), ass. (0->86)
t88 = cos(qJ(2));
t106 = t88 * qJD(1);
t84 = sin(qJ(2));
t107 = qJD(1) * t84;
t80 = sin(pkin(9));
t81 = cos(pkin(9));
t66 = t81 * qJD(2) - t80 * t107;
t67 = t80 * qJD(2) + t81 * t107;
t51 = -t66 * mrSges(4,1) + t67 * mrSges(4,2);
t52 = mrSges(4,2) * t106 + t66 * mrSges(4,3);
t105 = qJD(1) * qJD(2);
t101 = t88 * t105;
t71 = t84 * qJDD(1) + t101;
t55 = t80 * qJDD(2) + t81 * t71;
t77 = t84 * t105;
t72 = t88 * qJDD(1) - t77;
t85 = sin(qJ(1));
t89 = cos(qJ(1));
t103 = t85 * g(1) - t89 * g(2);
t91 = qJD(1) ^ 2;
t62 = -qJDD(1) * pkin(1) - t91 * pkin(6) - t103;
t40 = (-t71 - t101) * qJ(3) + (-t72 + t77) * pkin(2) + t62;
t98 = -t89 * g(1) - t85 * g(2);
t63 = -t91 * pkin(1) + qJDD(1) * pkin(6) + t98;
t102 = -t84 * g(3) + t88 * t63;
t69 = (-pkin(2) * t88 - qJ(3) * t84) * qJD(1);
t90 = qJD(2) ^ 2;
t43 = -t90 * pkin(2) + qJDD(2) * qJ(3) + t69 * t106 + t102;
t99 = -0.2e1 * qJD(3) * t67 + t81 * t40 - t80 * t43;
t22 = (-t66 * t106 - t55) * pkin(7) + (t66 * t67 - t72) * pkin(3) + t99;
t104 = 0.2e1 * qJD(3) * t66 + t80 * t40 + t81 * t43;
t54 = t81 * qJDD(2) - t80 * t71;
t56 = -pkin(3) * t106 - t67 * pkin(7);
t65 = t66 ^ 2;
t24 = -t65 * pkin(3) + t54 * pkin(7) + t56 * t106 + t104;
t83 = sin(qJ(4));
t87 = cos(qJ(4));
t100 = t87 * t22 - t83 * t24;
t49 = t87 * t66 - t83 * t67;
t33 = t49 * qJD(4) + t83 * t54 + t87 * t55;
t50 = t83 * t66 + t87 * t67;
t68 = qJDD(4) - t72;
t76 = qJD(4) - t106;
t13 = (t49 * t76 - t33) * pkin(8) + (t49 * t50 + t68) * pkin(4) + t100;
t109 = t83 * t22 + t87 * t24;
t32 = -t50 * qJD(4) + t87 * t54 - t83 * t55;
t46 = t76 * pkin(4) - t50 * pkin(8);
t48 = t49 ^ 2;
t14 = -t48 * pkin(4) + t32 * pkin(8) - t76 * t46 + t109;
t82 = sin(qJ(5));
t86 = cos(qJ(5));
t35 = t86 * t49 - t82 * t50;
t19 = t35 * qJD(5) + t82 * t32 + t86 * t33;
t36 = t82 * t49 + t86 * t50;
t26 = -t35 * mrSges(6,1) + t36 * mrSges(6,2);
t75 = qJD(5) + t76;
t29 = -t75 * mrSges(6,2) + t35 * mrSges(6,3);
t64 = qJDD(5) + t68;
t11 = m(6) * (t86 * t13 - t82 * t14) - t19 * mrSges(6,3) + t64 * mrSges(6,1) - t36 * t26 + t75 * t29;
t18 = -t36 * qJD(5) + t86 * t32 - t82 * t33;
t30 = t75 * mrSges(6,1) - t36 * mrSges(6,3);
t12 = m(6) * (t82 * t13 + t86 * t14) + t18 * mrSges(6,3) - t64 * mrSges(6,2) + t35 * t26 - t75 * t30;
t37 = -t49 * mrSges(5,1) + t50 * mrSges(5,2);
t44 = -t76 * mrSges(5,2) + t49 * mrSges(5,3);
t7 = m(5) * t100 + t68 * mrSges(5,1) - t33 * mrSges(5,3) + t86 * t11 + t82 * t12 - t50 * t37 + t76 * t44;
t45 = t76 * mrSges(5,1) - t50 * mrSges(5,3);
t8 = m(5) * t109 - t68 * mrSges(5,2) + t32 * mrSges(5,3) - t82 * t11 + t86 * t12 + t49 * t37 - t76 * t45;
t5 = m(4) * t99 - t72 * mrSges(4,1) - t55 * mrSges(4,3) - t52 * t106 - t67 * t51 + t87 * t7 + t83 * t8;
t53 = -mrSges(4,1) * t106 - t67 * mrSges(4,3);
t6 = m(4) * t104 + t72 * mrSges(4,2) + t54 * mrSges(4,3) + t53 * t106 + t66 * t51 - t83 * t7 + t87 * t8;
t73 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t107;
t74 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t106;
t111 = m(3) * t62 - t72 * mrSges(3,1) + t71 * mrSges(3,2) + t81 * t5 + t80 * t6 + (t84 * t73 - t88 * t74) * qJD(1);
t108 = -t88 * g(3) - t84 * t63;
t70 = (-mrSges(3,1) * t88 + mrSges(3,2) * t84) * qJD(1);
t42 = -qJDD(2) * pkin(2) - t90 * qJ(3) + t69 * t107 + qJDD(3) - t108;
t94 = -t54 * pkin(3) - t65 * pkin(7) + t67 * t56 + t42;
t96 = t18 * mrSges(6,1) + t35 * t29 - m(6) * (-t32 * pkin(4) - t48 * pkin(8) + t50 * t46 + t94) - t19 * mrSges(6,2) - t36 * t30;
t93 = -m(5) * t94 + t32 * mrSges(5,1) - t33 * mrSges(5,2) + t49 * t44 - t50 * t45 + t96;
t92 = m(4) * t42 - t54 * mrSges(4,1) + t55 * mrSges(4,2) - t66 * t52 + t67 * t53 - t93;
t10 = m(3) * t108 + qJDD(2) * mrSges(3,1) - t71 * mrSges(3,3) + qJD(2) * t74 - t70 * t107 - t92;
t4 = m(3) * t102 - qJDD(2) * mrSges(3,2) + t72 * mrSges(3,3) - qJD(2) * t73 + t70 * t106 - t80 * t5 + t81 * t6;
t110 = t88 * t10 + t84 * t4;
t2 = m(2) * t103 + qJDD(1) * mrSges(2,1) - t91 * mrSges(2,2) - t111;
t1 = m(2) * t98 - t91 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t84 * t10 + t88 * t4;
t3 = [-m(1) * g(1) + t89 * t1 - t85 * t2, t1, t4, t6, t8, t12; -m(1) * g(2) + t85 * t1 + t89 * t2, t2, t10, t5, t7, t11; (-m(1) - m(2)) * g(3) + t110, -m(2) * g(3) + t110, t111, t92, -t93, -t96;];
f_new = t3;
