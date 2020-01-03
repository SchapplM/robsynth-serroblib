% Calculate vector of cutting forces with Newton-Euler
% S5RPPRP5
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
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
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPPRP5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP5_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP5_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP5_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP5_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP5_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:22
% EndTime: 2019-12-31 17:53:23
% DurationCPUTime: 0.52s
% Computational Cost: add. (3189->127), mult. (7537->163), div. (0->0), fcn. (4413->6), ass. (0->69)
t115 = mrSges(4,2) + mrSges(3,3);
t69 = sin(pkin(7));
t100 = qJD(1) * t69;
t101 = qJ(3) * t69;
t72 = sin(qJ(1));
t74 = cos(qJ(1));
t103 = t72 * g(1) - t74 * g(2);
t76 = qJD(1) ^ 2;
t47 = -qJDD(1) * pkin(1) - t76 * qJ(2) + qJDD(2) - t103;
t70 = cos(pkin(7));
t96 = qJDD(1) * t70;
t114 = -pkin(2) * t96 - 0.2e1 * qJD(3) * t100 - qJDD(1) * t101 + t47;
t65 = t70 ^ 2;
t102 = t69 ^ 2 + t65;
t93 = t102 * t76;
t113 = (mrSges(3,2) - mrSges(4,3)) * t69 - (mrSges(3,1) + mrSges(4,1)) * t70;
t89 = -t74 * g(1) - t72 * g(2);
t50 = -t76 * pkin(1) + qJDD(1) * qJ(2) + t89;
t104 = -t70 * g(3) - t69 * t50;
t88 = -t70 * mrSges(4,1) - t69 * mrSges(4,3);
t52 = t88 * qJD(1);
t53 = (-t70 * mrSges(3,1) + t69 * mrSges(3,2)) * qJD(1);
t51 = (-pkin(2) * t70 - t101) * qJD(1);
t94 = 0.2e1 * qJD(1) * qJD(2);
t20 = t51 * t100 + t69 * t94 + qJDD(3) - t104;
t71 = sin(qJ(4));
t73 = cos(qJ(4));
t85 = t69 * t71 + t70 * t73;
t48 = t85 * qJD(1);
t86 = t69 * t73 - t70 * t71;
t49 = t86 * qJD(1);
t28 = t48 * mrSges(6,1) - t49 * mrSges(6,3);
t105 = -t48 * mrSges(5,1) - t49 * mrSges(5,2) - t28;
t110 = pkin(3) * t76;
t15 = (-pkin(6) * qJDD(1) - t70 * t110) * t69 + t20;
t91 = -t69 * g(3) + (t50 + t94) * t70;
t99 = qJD(1) * t70;
t84 = t51 * t99 + t91;
t17 = -pkin(6) * t96 - t65 * t110 + t84;
t106 = t71 * t15 + t73 * t17;
t107 = -mrSges(5,3) - mrSges(6,2);
t97 = t49 * qJD(4);
t35 = t85 * qJDD(1) + t97;
t39 = qJD(4) * mrSges(5,1) - t49 * mrSges(5,3);
t27 = t48 * pkin(4) - t49 * qJ(5);
t40 = -qJD(4) * mrSges(6,1) + t49 * mrSges(6,2);
t75 = qJD(4) ^ 2;
t95 = m(6) * (-t75 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) - t48 * t27 + t106) + qJD(4) * t40 + qJDD(4) * mrSges(6,3);
t8 = m(5) * t106 - qJDD(4) * mrSges(5,2) - qJD(4) * t39 + t105 * t48 + t107 * t35 + t95;
t87 = t73 * t15 - t71 * t17;
t111 = m(6) * (-qJDD(4) * pkin(4) - t75 * qJ(5) + t49 * t27 + qJDD(5) - t87);
t98 = t48 * qJD(4);
t36 = t86 * qJDD(1) - t98;
t38 = -qJD(4) * mrSges(5,2) - t48 * mrSges(5,3);
t41 = -t48 * mrSges(6,2) + qJD(4) * mrSges(6,3);
t9 = m(5) * t87 - t111 + t105 * t49 + t107 * t36 + (mrSges(5,1) + mrSges(6,1)) * qJDD(4) + (t38 + t41) * qJD(4);
t81 = m(4) * t20 + t71 * t8 + t73 * t9;
t4 = m(3) * t104 + (-t115 * qJDD(1) + (-0.2e1 * m(3) * qJD(2) - t52 - t53) * qJD(1)) * t69 - t81;
t83 = m(4) * t84 + mrSges(4,2) * t96 + t52 * t99 - t71 * t9 + t73 * t8;
t5 = m(3) * t91 + (qJDD(1) * mrSges(3,3) + qJD(1) * t53) * t70 + t83;
t112 = t70 * t4 + t69 * t5;
t77 = pkin(3) * t96 - pkin(6) * t93 - t114;
t90 = m(6) * (-0.2e1 * qJD(5) * t49 + (-t36 + t98) * qJ(5) + (t35 + t97) * pkin(4) + t77) - t36 * mrSges(6,3) + t35 * mrSges(6,1) - t49 * t40 + t48 * t41;
t80 = m(5) * t77 + t35 * mrSges(5,1) + t36 * mrSges(5,2) + t48 * t38 + t49 * t39 + t90;
t79 = m(4) * t114 - t80;
t78 = m(3) * t47 + t79;
t6 = -t78 + (mrSges(2,1) - t113) * qJDD(1) + (t115 * t102 - mrSges(2,2)) * t76 + m(2) * t103;
t1 = m(2) * t89 - t76 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t69 * t4 + t70 * t5;
t2 = [-m(1) * g(1) + t74 * t1 - t72 * t6, t1, t5, t83, t8, -t35 * mrSges(6,2) - t48 * t28 + t95; -m(1) * g(2) + t72 * t1 + t74 * t6, t6, t4, -mrSges(4,2) * t93 + t88 * qJDD(1) + t79, t9, t90; (-m(1) - m(2)) * g(3) + t112, -m(2) * g(3) + t112, t113 * qJDD(1) - t115 * t93 + t78, (qJDD(1) * mrSges(4,2) + qJD(1) * t52) * t69 + t81, t80, -qJDD(4) * mrSges(6,1) + t36 * mrSges(6,2) - qJD(4) * t41 + t49 * t28 + t111;];
f_new = t2;
