% Calculate vector of cutting forces with Newton-Euler
% S5RRRRR9
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRRR9_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR9_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR9_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR9_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR9_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR9_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR9_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR9_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR9_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:27:49
% EndTime: 2019-12-31 22:27:56
% DurationCPUTime: 2.53s
% Computational Cost: add. (32805->165), mult. (66550->216), div. (0->0), fcn. (46468->10), ass. (0->88)
t83 = sin(qJ(2));
t106 = qJD(1) * t83;
t82 = sin(qJ(3));
t87 = cos(qJ(3));
t66 = t87 * qJD(2) - t82 * t106;
t104 = qJD(1) * qJD(2);
t88 = cos(qJ(2));
t101 = t88 * t104;
t70 = t83 * qJDD(1) + t101;
t49 = t66 * qJD(3) + t82 * qJDD(2) + t87 * t70;
t67 = t82 * qJD(2) + t87 * t106;
t53 = -t66 * mrSges(4,1) + t67 * mrSges(4,2);
t105 = t88 * qJD(1);
t76 = qJD(3) - t105;
t54 = -t76 * mrSges(4,2) + t66 * mrSges(4,3);
t77 = t83 * t104;
t71 = t88 * qJDD(1) - t77;
t65 = qJDD(3) - t71;
t84 = sin(qJ(1));
t89 = cos(qJ(1));
t103 = t84 * g(1) - t89 * g(2);
t91 = qJD(1) ^ 2;
t61 = -qJDD(1) * pkin(1) - t91 * pkin(6) - t103;
t40 = (-t70 - t101) * pkin(7) + (-t71 + t77) * pkin(2) + t61;
t98 = -t89 * g(1) - t84 * g(2);
t62 = -t91 * pkin(1) + qJDD(1) * pkin(6) + t98;
t102 = -t83 * g(3) + t88 * t62;
t69 = (-pkin(2) * t88 - pkin(7) * t83) * qJD(1);
t90 = qJD(2) ^ 2;
t43 = -t90 * pkin(2) + qJDD(2) * pkin(7) + t69 * t105 + t102;
t99 = t87 * t40 - t82 * t43;
t22 = (t66 * t76 - t49) * pkin(8) + (t66 * t67 + t65) * pkin(3) + t99;
t108 = t82 * t40 + t87 * t43;
t48 = -t67 * qJD(3) + t87 * qJDD(2) - t82 * t70;
t56 = t76 * pkin(3) - t67 * pkin(8);
t64 = t66 ^ 2;
t24 = -t64 * pkin(3) + t48 * pkin(8) - t76 * t56 + t108;
t81 = sin(qJ(4));
t86 = cos(qJ(4));
t100 = t86 * t22 - t81 * t24;
t51 = t86 * t66 - t81 * t67;
t31 = t51 * qJD(4) + t81 * t48 + t86 * t49;
t52 = t81 * t66 + t86 * t67;
t63 = qJDD(4) + t65;
t75 = qJD(4) + t76;
t13 = (t51 * t75 - t31) * pkin(9) + (t51 * t52 + t63) * pkin(4) + t100;
t109 = t81 * t22 + t86 * t24;
t30 = -t52 * qJD(4) + t86 * t48 - t81 * t49;
t46 = t75 * pkin(4) - t52 * pkin(9);
t50 = t51 ^ 2;
t14 = -t50 * pkin(4) + t30 * pkin(9) - t75 * t46 + t109;
t80 = sin(qJ(5));
t85 = cos(qJ(5));
t35 = t85 * t51 - t80 * t52;
t19 = t35 * qJD(5) + t80 * t30 + t85 * t31;
t36 = t80 * t51 + t85 * t52;
t26 = -t35 * mrSges(6,1) + t36 * mrSges(6,2);
t72 = qJD(5) + t75;
t32 = -t72 * mrSges(6,2) + t35 * mrSges(6,3);
t60 = qJDD(5) + t63;
t11 = m(6) * (t85 * t13 - t80 * t14) - t19 * mrSges(6,3) + t60 * mrSges(6,1) - t36 * t26 + t72 * t32;
t18 = -t36 * qJD(5) + t85 * t30 - t80 * t31;
t33 = t72 * mrSges(6,1) - t36 * mrSges(6,3);
t12 = m(6) * (t80 * t13 + t85 * t14) + t18 * mrSges(6,3) - t60 * mrSges(6,2) + t35 * t26 - t72 * t33;
t37 = -t51 * mrSges(5,1) + t52 * mrSges(5,2);
t44 = -t75 * mrSges(5,2) + t51 * mrSges(5,3);
t7 = m(5) * t100 + t63 * mrSges(5,1) - t31 * mrSges(5,3) + t85 * t11 + t80 * t12 - t52 * t37 + t75 * t44;
t45 = t75 * mrSges(5,1) - t52 * mrSges(5,3);
t8 = m(5) * t109 - t63 * mrSges(5,2) + t30 * mrSges(5,3) - t80 * t11 + t85 * t12 + t51 * t37 - t75 * t45;
t5 = m(4) * t99 + t65 * mrSges(4,1) - t49 * mrSges(4,3) - t67 * t53 + t76 * t54 + t86 * t7 + t81 * t8;
t55 = t76 * mrSges(4,1) - t67 * mrSges(4,3);
t6 = m(4) * t108 - t65 * mrSges(4,2) + t48 * mrSges(4,3) + t66 * t53 - t76 * t55 - t81 * t7 + t86 * t8;
t73 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t106;
t74 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t105;
t111 = m(3) * t61 - t71 * mrSges(3,1) + t70 * mrSges(3,2) + t87 * t5 + t82 * t6 + (t83 * t73 - t88 * t74) * qJD(1);
t107 = -t88 * g(3) - t83 * t62;
t68 = (-mrSges(3,1) * t88 + mrSges(3,2) * t83) * qJD(1);
t42 = -qJDD(2) * pkin(2) - t90 * pkin(7) + t69 * t106 - t107;
t94 = -t48 * pkin(3) - t64 * pkin(8) + t67 * t56 + t42;
t96 = t18 * mrSges(6,1) + t35 * t32 - m(6) * (-t30 * pkin(4) - t50 * pkin(9) + t52 * t46 + t94) - t19 * mrSges(6,2) - t36 * t33;
t93 = -m(5) * t94 + t30 * mrSges(5,1) - t31 * mrSges(5,2) + t51 * t44 - t52 * t45 + t96;
t92 = m(4) * t42 - t48 * mrSges(4,1) + t49 * mrSges(4,2) - t66 * t54 + t67 * t55 - t93;
t10 = m(3) * t107 + qJDD(2) * mrSges(3,1) - t70 * mrSges(3,3) + qJD(2) * t74 - t68 * t106 - t92;
t4 = m(3) * t102 - qJDD(2) * mrSges(3,2) + t71 * mrSges(3,3) - qJD(2) * t73 + t68 * t105 - t82 * t5 + t87 * t6;
t110 = t88 * t10 + t83 * t4;
t2 = m(2) * t103 + qJDD(1) * mrSges(2,1) - t91 * mrSges(2,2) - t111;
t1 = m(2) * t98 - t91 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t83 * t10 + t88 * t4;
t3 = [-m(1) * g(1) + t89 * t1 - t84 * t2, t1, t4, t6, t8, t12; -m(1) * g(2) + t84 * t1 + t89 * t2, t2, t10, t5, t7, t11; (-m(1) - m(2)) * g(3) + t110, -m(2) * g(3) + t110, t111, t92, -t93, -t96;];
f_new = t3;
