% Calculate vector of cutting forces with Newton-Euler
% S5RPPPR6
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
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
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPPPR6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR6_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR6_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR6_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR6_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR6_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:47:39
% EndTime: 2019-12-31 17:47:40
% DurationCPUTime: 0.77s
% Computational Cost: add. (4448->140), mult. (11452->194), div. (0->0), fcn. (6531->8), ass. (0->86)
t116 = 2 * qJD(1) * qJD(2);
t62 = sin(pkin(7));
t100 = qJ(3) * t62;
t64 = cos(pkin(7));
t80 = -pkin(2) * t64 - t100;
t45 = t80 * qJD(1);
t69 = qJD(1) ^ 2;
t66 = sin(qJ(1));
t68 = cos(qJ(1));
t82 = -t68 * g(1) - t66 * g(2);
t44 = -(t69 * pkin(1)) + qJDD(1) * qJ(2) + t82;
t84 = -t62 * g(3) + (t116 + t44) * t64;
t97 = qJD(1) * t64;
t115 = -t45 * t97 - t84;
t114 = mrSges(4,1) + mrSges(3,3);
t104 = t64 * mrSges(3,1);
t113 = (mrSges(3,2) - mrSges(4,3)) * t62 - t104;
t112 = -2 * qJD(4);
t102 = -t64 * g(3) - t62 * t44;
t46 = (mrSges(4,2) * t64 - t62 * mrSges(4,3)) * qJD(1);
t47 = (t62 * mrSges(3,2) - t104) * qJD(1);
t98 = qJD(1) * t62;
t23 = t62 * t116 + t45 * t98 + qJDD(3) - t102;
t61 = sin(pkin(8));
t63 = cos(pkin(8));
t105 = t63 * t64;
t59 = t62 ^ 2;
t109 = t59 * t69;
t38 = (pkin(4) * t63 + pkin(6) * t61) * t97;
t92 = t63 * t97;
t99 = qJ(4) * t69;
t18 = (pkin(3) * qJDD(1) - t64 * t99) * t62 + t23;
t60 = t64 ^ 2;
t101 = t59 + t60;
t52 = -0.2e1 * qJD(3) * t98;
t90 = t66 * g(1) - t68 * g(2);
t81 = qJDD(2) - t90;
t22 = t52 + (-t101 * pkin(3) - qJ(2)) * t69 + (-t100 - pkin(1) + (-pkin(2) - qJ(4)) * t64) * qJDD(1) + t81;
t93 = t92 * t112 + t61 * t18 + t63 * t22;
t96 = qJDD(1) * t62;
t15 = -pkin(4) * t109 + pkin(6) * t96 - t38 * t92 + t93;
t107 = t62 * t69;
t95 = qJDD(1) * t64;
t71 = pkin(3) * t95 - t60 * t99 + qJDD(4) - t115;
t16 = ((qJDD(1) * t61 + t63 * t107) * pkin(6) + (qJDD(1) * t63 - t61 * t107) * pkin(4)) * t64 + t71;
t108 = t61 * t64;
t65 = sin(qJ(5));
t67 = cos(qJ(5));
t76 = t65 * t108 + t62 * t67;
t35 = t76 * qJD(1);
t77 = -t67 * t108 + t62 * t65;
t36 = t77 * qJD(1);
t25 = -t35 * mrSges(6,1) + t36 * mrSges(6,2);
t27 = t35 * qJD(5) + t77 * qJDD(1);
t49 = qJD(5) + t92;
t30 = -t49 * mrSges(6,2) + t35 * mrSges(6,3);
t91 = t63 * t95;
t48 = qJDD(5) + t91;
t12 = m(6) * (-t65 * t15 + t67 * t16) - t27 * mrSges(6,3) + t48 * mrSges(6,1) - t36 * t25 + t49 * t30;
t26 = -t36 * qJD(5) + t76 * qJDD(1);
t31 = t49 * mrSges(6,1) - t36 * mrSges(6,3);
t13 = m(6) * (t67 * t15 + t65 * t16) + t26 * mrSges(6,3) - t48 * mrSges(6,2) + t35 * t25 - t49 * t31;
t110 = mrSges(5,2) * t61;
t37 = (mrSges(5,1) * t63 - t110) * t97;
t42 = (t62 * mrSges(5,1) + mrSges(5,3) * t108) * qJD(1);
t78 = -t62 * mrSges(5,2) - mrSges(5,3) * t105;
t8 = m(5) * t93 + t67 * t13 - t65 * t12 + t78 * qJDD(1) + (-t37 * t105 - t62 * t42) * qJD(1);
t106 = t63 * t18;
t43 = t78 * qJD(1);
t70 = m(6) * (-pkin(4) * t96 - pkin(6) * t109 - t106 + (t22 + (t112 - t38) * t97) * t61) - t26 * mrSges(6,1) + t27 * mrSges(6,2) - t35 * t30 + t36 * t31;
t9 = m(5) * t106 + (qJDD(1) * mrSges(5,1) + qJD(1) * t43) * t62 + (-m(5) * t22 + (qJDD(1) * mrSges(5,3) + (0.2e1 * m(5) * qJD(4) + t37) * qJD(1)) * t64) * t61 - t70;
t75 = m(4) * t23 + t61 * t8 + t63 * t9;
t3 = m(3) * t102 + (-t114 * qJDD(1) + (-0.2e1 * m(3) * qJD(2) - t46 - t47) * qJD(1)) * t62 - t75;
t83 = m(5) * t71 + mrSges(5,1) * t91 + t67 * t12 + t65 * t13 + t43 * t92;
t74 = m(4) * t115 - t83;
t85 = t42 * t61 - t46;
t86 = -mrSges(4,1) + t110;
t6 = m(3) * t84 + ((mrSges(3,3) - t86) * qJDD(1) + (t47 - t85) * qJD(1)) * t64 - t74;
t111 = t64 * t3 + t62 * t6;
t88 = t101 * mrSges(4,1);
t73 = -t69 * qJ(2) + t81;
t79 = -t61 * t9 + t63 * t8 + m(4) * (t52 + (-pkin(1) + t80) * qJDD(1) + t73) + mrSges(4,2) * t95;
t72 = m(3) * (-qJDD(1) * pkin(1) + t73) + t79;
t4 = m(2) * t90 + (mrSges(2,1) - t113) * qJDD(1) + (t114 * t101 - mrSges(2,2)) * t69 - t72;
t1 = m(2) * t82 - t69 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t62 * t3 + t64 * t6;
t2 = [-m(1) * g(1) + t68 * t1 - t66 * t4, t1, t6, -mrSges(4,3) * t96 - t69 * t88 + t79, t8, t13; -m(1) * g(2) + t66 * t1 + t68 * t4, t4, t3, (t85 * qJD(1) + t86 * qJDD(1)) * t64 + t74, t9, t12; (-m(1) - m(2)) * g(3) + t111, -m(2) * g(3) + t111, t113 * qJDD(1) + (-t101 * mrSges(3,3) - t88) * t69 + t72, (qJDD(1) * mrSges(4,1) + qJD(1) * t46) * t62 + t75, (-qJDD(1) * mrSges(5,2) - qJD(1) * t42) * t108 + t83, t70;];
f_new = t2;
