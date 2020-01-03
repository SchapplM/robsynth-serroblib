% Calculate vector of cutting forces with Newton-Euler
% S5RPRRR9
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRRR9_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR9_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR9_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR9_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR9_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR9_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR9_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR9_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR9_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:07:23
% EndTime: 2019-12-31 19:07:27
% DurationCPUTime: 1.95s
% Computational Cost: add. (22179->154), mult. (53429->201), div. (0->0), fcn. (40273->10), ass. (0->85)
t84 = qJD(1) ^ 2;
t75 = cos(pkin(9));
t72 = t75 ^ 2;
t74 = sin(pkin(9));
t105 = t74 ^ 2 + t72;
t110 = mrSges(3,3) * t105;
t102 = qJD(1) * qJD(2);
t100 = -t75 * g(3) - 0.2e1 * t74 * t102;
t78 = sin(qJ(3));
t82 = cos(qJ(3));
t92 = -t74 * t78 + t75 * t82;
t62 = t92 * qJD(1);
t93 = t74 * t82 + t75 * t78;
t63 = t93 * qJD(1);
t77 = sin(qJ(4));
t81 = cos(qJ(4));
t46 = t81 * t62 - t77 * t63;
t54 = -t63 * qJD(3) + qJDD(1) * t92;
t103 = t62 * qJD(3);
t55 = qJDD(1) * t93 + t103;
t29 = t46 * qJD(4) + t77 * t54 + t81 * t55;
t47 = t77 * t62 + t81 * t63;
t34 = -t46 * mrSges(5,1) + t47 * mrSges(5,2);
t73 = qJD(3) + qJD(4);
t40 = -t73 * mrSges(5,2) + t46 * mrSges(5,3);
t70 = qJDD(3) + qJDD(4);
t76 = sin(qJ(5));
t80 = cos(qJ(5));
t37 = t80 * t47 + t76 * t73;
t20 = -t37 * qJD(5) - t76 * t29 + t80 * t70;
t36 = -t76 * t47 + t80 * t73;
t21 = t36 * qJD(5) + t80 * t29 + t76 * t70;
t45 = qJD(5) - t46;
t31 = -t45 * mrSges(6,2) + t36 * mrSges(6,3);
t32 = t45 * mrSges(6,1) - t37 * mrSges(6,3);
t35 = -t46 * pkin(4) - t47 * pkin(8);
t69 = t73 ^ 2;
t104 = pkin(6) * qJDD(1);
t108 = pkin(2) * t84;
t79 = sin(qJ(1));
t83 = cos(qJ(1));
t97 = -t83 * g(1) - t79 * g(2);
t64 = -t84 * pkin(1) + qJDD(1) * qJ(2) + t97;
t43 = (t108 * t75 - t104 - t64) * t74 + t100;
t98 = -t74 * g(3) + (0.2e1 * t102 + t64) * t75;
t44 = t104 * t75 - t108 * t72 + t98;
t99 = t82 * t43 - t78 * t44;
t19 = (-t55 + t103) * pkin(7) + (t62 * t63 + qJDD(3)) * pkin(3) + t99;
t106 = t78 * t43 + t82 * t44;
t58 = qJD(3) * pkin(3) - t63 * pkin(7);
t61 = t62 ^ 2;
t23 = -t61 * pkin(3) + t54 * pkin(7) - qJD(3) * t58 + t106;
t94 = t81 * t19 - t77 * t23;
t89 = m(6) * (-t70 * pkin(4) - t69 * pkin(8) + t47 * t35 - t94) - t20 * mrSges(6,1) + t21 * mrSges(6,2) - t36 * t31 + t37 * t32;
t10 = m(5) * t94 + t70 * mrSges(5,1) - t29 * mrSges(5,3) - t47 * t34 + t73 * t40 - t89;
t51 = -t62 * mrSges(4,1) + t63 * mrSges(4,2);
t56 = -qJD(3) * mrSges(4,2) + t62 * mrSges(4,3);
t107 = t77 * t19 + t81 * t23;
t16 = -t69 * pkin(4) + t70 * pkin(8) + t46 * t35 + t107;
t28 = -t47 * qJD(4) + t81 * t54 - t77 * t55;
t101 = t79 * g(1) - t83 * g(2);
t96 = qJDD(2) - t101;
t87 = (-pkin(2) * t75 - pkin(1)) * qJDD(1) + (-pkin(6) * t105 - qJ(2)) * t84 + t96;
t85 = -t54 * pkin(3) - t61 * pkin(7) + t63 * t58 + t87;
t17 = (-t46 * t73 - t29) * pkin(8) + (t47 * t73 - t28) * pkin(4) + t85;
t27 = qJDD(5) - t28;
t30 = -t36 * mrSges(6,1) + t37 * mrSges(6,2);
t13 = m(6) * (-t76 * t16 + t80 * t17) - t21 * mrSges(6,3) + t27 * mrSges(6,1) - t37 * t30 + t45 * t31;
t14 = m(6) * (t80 * t16 + t76 * t17) + t20 * mrSges(6,3) - t27 * mrSges(6,2) + t36 * t30 - t45 * t32;
t41 = t73 * mrSges(5,1) - t47 * mrSges(5,3);
t9 = m(5) * t107 - t70 * mrSges(5,2) + t28 * mrSges(5,3) - t76 * t13 + t80 * t14 + t46 * t34 - t73 * t41;
t6 = m(4) * t99 + qJDD(3) * mrSges(4,1) - t55 * mrSges(4,3) + qJD(3) * t56 + t81 * t10 - t63 * t51 + t77 * t9;
t57 = qJD(3) * mrSges(4,1) - t63 * mrSges(4,3);
t7 = m(4) * t106 - qJDD(3) * mrSges(4,2) + t54 * mrSges(4,3) - qJD(3) * t57 - t77 * t10 + t62 * t51 + t81 * t9;
t95 = -t75 * mrSges(3,1) + t74 * mrSges(3,2);
t91 = qJDD(1) * mrSges(3,3) + t84 * t95;
t4 = m(3) * t100 + t78 * t7 + t82 * t6 + (-m(3) * t64 - t91) * t74;
t5 = m(3) * t98 - t78 * t6 + t82 * t7 + t75 * t91;
t109 = t75 * t4 + t74 * t5;
t90 = -m(5) * t85 + t28 * mrSges(5,1) - t29 * mrSges(5,2) - t80 * t13 - t76 * t14 + t46 * t40 - t47 * t41;
t88 = -m(4) * t87 + t54 * mrSges(4,1) - t55 * mrSges(4,2) + t62 * t56 - t63 * t57 + t90;
t86 = m(3) * (-qJDD(1) * pkin(1) - t84 * qJ(2) + t96) - t88;
t8 = -t86 + (-mrSges(2,2) + t110) * t84 + (mrSges(2,1) - t95) * qJDD(1) + m(2) * t101;
t1 = m(2) * t97 - t84 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t74 * t4 + t75 * t5;
t2 = [-m(1) * g(1) + t83 * t1 - t79 * t8, t1, t5, t7, t9, t14; -m(1) * g(2) + t79 * t1 + t83 * t8, t8, t4, t6, t10, t13; (-m(1) - m(2)) * g(3) + t109, -m(2) * g(3) + t109, qJDD(1) * t95 - t84 * t110 + t86, -t88, -t90, t89;];
f_new = t2;
