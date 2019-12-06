% Calculate vector of cutting forces with Newton-Euler
% S5RRRRR3
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
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
% Datum: 2019-12-05 18:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRRR3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR3_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR3_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_invdynf_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR3_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR3_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR3_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:55:10
% EndTime: 2019-12-05 18:55:15
% DurationCPUTime: 1.78s
% Computational Cost: add. (21353->152), mult. (43366->203), div. (0->0), fcn. (32558->10), ass. (0->81)
t101 = qJD(1) * qJD(2);
t82 = sin(qJ(2));
t87 = cos(qJ(2));
t67 = t82 * qJDD(1) + t87 * t101;
t99 = t82 * t101;
t68 = t87 * qJDD(1) - t99;
t103 = qJD(1) * t82;
t69 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t103;
t102 = qJD(1) * t87;
t70 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t102;
t81 = sin(qJ(3));
t86 = cos(qJ(3));
t63 = (t81 * t87 + t82 * t86) * qJD(1);
t41 = -t63 * qJD(3) - t81 * t67 + t86 * t68;
t40 = qJDD(4) - t41;
t78 = qJD(2) + qJD(3);
t80 = sin(qJ(4));
t85 = cos(qJ(4));
t54 = -t80 * t63 + t85 * t78;
t55 = t85 * t63 + t80 * t78;
t62 = t86 * t102 - t81 * t103;
t42 = t62 * qJD(3) + t86 * t67 + t81 * t68;
t83 = sin(qJ(1));
t88 = cos(qJ(1));
t71 = t83 * g(1) - t88 * g(2);
t94 = -t71 + (-t68 + t99) * pkin(1);
t24 = (-t62 * t78 - t42) * pkin(5) + (t63 * t78 - t41) * pkin(2) + t94;
t89 = qJD(1) ^ 2;
t72 = -t88 * g(1) - t83 * g(2);
t96 = -t87 * g(3) - t82 * t72;
t53 = (t82 * t87 * t89 + qJDD(2)) * pkin(1) + t96;
t100 = -t82 * g(3) + t87 * t72;
t56 = (-t87 ^ 2 * t89 - qJD(2) ^ 2) * pkin(1) + t100;
t104 = t81 * t53 + t86 * t56;
t47 = -t62 * pkin(2) - t63 * pkin(5);
t76 = t78 ^ 2;
t77 = qJDD(2) + qJDD(3);
t29 = -t76 * pkin(2) + t77 * pkin(5) + t62 * t47 + t104;
t98 = t85 * t24 - t80 * t29;
t15 = (t54 * t55 + t40) * pkin(3) + t98;
t105 = t80 * t24 + t85 * t29;
t61 = qJD(4) - t62;
t16 = (-t54 ^ 2 - t61 ^ 2) * pkin(3) + t105;
t30 = -t55 * qJD(4) - t80 * t42 + t85 * t77;
t31 = t54 * qJD(4) + t85 * t42 + t80 * t77;
t79 = sin(qJ(5));
t84 = cos(qJ(5));
t34 = t84 * t54 - t79 * t55;
t21 = t34 * qJD(5) + t79 * t30 + t84 * t31;
t35 = t79 * t54 + t84 * t55;
t27 = -t34 * mrSges(6,1) + t35 * mrSges(6,2);
t59 = qJD(5) + t61;
t32 = -t59 * mrSges(6,2) + t34 * mrSges(6,3);
t37 = qJDD(5) + t40;
t13 = m(6) * (t84 * t15 - t79 * t16) - t21 * mrSges(6,3) + t37 * mrSges(6,1) - t35 * t27 + t59 * t32;
t20 = -t35 * qJD(5) + t84 * t30 - t79 * t31;
t33 = t59 * mrSges(6,1) - t35 * mrSges(6,3);
t14 = m(6) * (t79 * t15 + t84 * t16) + t20 * mrSges(6,3) - t37 * mrSges(6,2) + t34 * t27 - t59 * t33;
t36 = -t54 * mrSges(5,1) + t55 * mrSges(5,2);
t43 = -t61 * mrSges(5,2) + t54 * mrSges(5,3);
t10 = m(5) * t98 + t40 * mrSges(5,1) - t31 * mrSges(5,3) + t84 * t13 + t79 * t14 - t55 * t36 + t61 * t43;
t44 = t61 * mrSges(5,1) - t55 * mrSges(5,3);
t11 = m(5) * t105 - t40 * mrSges(5,2) + t30 * mrSges(5,3) - t79 * t13 + t84 * t14 + t54 * t36 - t61 * t44;
t57 = -t78 * mrSges(4,2) + t62 * mrSges(4,3);
t58 = t78 * mrSges(4,1) - t63 * mrSges(4,3);
t92 = -m(4) * t94 + t41 * mrSges(4,1) - t42 * mrSges(4,2) - t85 * t10 - t80 * t11 + t62 * t57 - t63 * t58;
t107 = t68 * mrSges(3,1) - t67 * mrSges(3,2) - (t69 * t82 - t70 * t87) * qJD(1) + t92;
t46 = -t62 * mrSges(4,1) + t63 * mrSges(4,2);
t97 = t86 * t53 - t81 * t56;
t28 = -t77 * pkin(2) - t76 * pkin(5) + t63 * t47 - t97;
t93 = t20 * mrSges(6,1) + t34 * t32 - m(6) * ((t55 * t61 - t30) * pkin(3) + t28) - t21 * mrSges(6,2) - t35 * t33;
t90 = m(5) * t28 - t30 * mrSges(5,1) + t31 * mrSges(5,2) - t54 * t43 + t55 * t44 - t93;
t12 = m(4) * t97 + t77 * mrSges(4,1) - t42 * mrSges(4,3) - t63 * t46 + t78 * t57 - t90;
t66 = (-mrSges(3,1) * t87 + mrSges(3,2) * t82) * qJD(1);
t7 = m(4) * t104 - t77 * mrSges(4,2) + t41 * mrSges(4,3) - t80 * t10 + t85 * t11 + t62 * t46 - t78 * t58;
t4 = m(3) * t96 + qJDD(2) * mrSges(3,1) - t67 * mrSges(3,3) + qJD(2) * t70 - t66 * t103 + t86 * t12 + t81 * t7;
t5 = m(3) * t100 - qJDD(2) * mrSges(3,2) + t68 * mrSges(3,3) - qJD(2) * t69 + t66 * t102 - t81 * t12 + t86 * t7;
t106 = t87 * t4 + t82 * t5;
t6 = qJDD(1) * mrSges(2,1) - t89 * mrSges(2,2) + (m(2) + m(3)) * t71 + t107;
t1 = m(2) * t72 - t89 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t82 * t4 + t87 * t5;
t2 = [-m(1) * g(1) + t88 * t1 - t83 * t6, t1, t5, t7, t11, t14; -m(1) * g(2) + t83 * t1 + t88 * t6, t6, t4, t12, t10, t13; (-m(1) - m(2)) * g(3) + t106, -m(2) * g(3) + t106, -m(3) * t71 - t107, -t92, t90, -t93;];
f_new = t2;
