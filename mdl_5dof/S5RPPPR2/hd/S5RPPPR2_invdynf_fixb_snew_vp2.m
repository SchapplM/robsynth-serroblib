% Calculate vector of cutting forces with Newton-Euler
% S5RPPPR2
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPPPR2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR2_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR2_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR2_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:22:25
% EndTime: 2020-01-03 11:22:28
% DurationCPUTime: 1.38s
% Computational Cost: add. (10432->149), mult. (29290->220), div. (0->0), fcn. (19772->10), ass. (0->93)
t118 = -2 * qJD(3);
t72 = sin(pkin(7));
t106 = qJD(1) * t72;
t77 = sin(qJ(1));
t79 = cos(qJ(1));
t108 = -t79 * g(2) - t77 * g(3);
t80 = qJD(1) ^ 2;
t85 = -t80 * qJ(2) + qJDD(2) - t108;
t75 = cos(pkin(7));
t90 = -pkin(2) * t75 - qJ(3) * t72;
t37 = (-pkin(1) + t90) * qJDD(1) + t85;
t117 = t106 * t118 + t37;
t68 = t72 ^ 2;
t69 = t75 ^ 2;
t116 = (t68 + t69) * mrSges(3,3);
t74 = cos(pkin(8));
t111 = t72 * t74;
t70 = sin(pkin(9));
t73 = cos(pkin(9));
t86 = t111 * t70 + t73 * t75;
t44 = t86 * qJD(1);
t42 = t86 * qJDD(1);
t115 = 2 * qJD(4);
t92 = -mrSges(3,1) * t75 + mrSges(3,2) * t72;
t56 = t92 * qJD(1);
t71 = sin(pkin(8));
t112 = t71 * t72;
t105 = qJD(1) * t75;
t98 = t74 * t106;
t45 = -t105 * t70 + t73 * t98;
t30 = mrSges(5,1) * t44 + mrSges(5,2) * t45;
t99 = t71 * t106;
t35 = -mrSges(5,2) * t99 - mrSges(5,3) * t44;
t43 = (t111 * t73 - t70 * t75) * qJDD(1);
t102 = t71 ^ 2 * t68 * t80;
t76 = sin(qJ(5));
t78 = cos(qJ(5));
t33 = t45 * t78 + t76 * t99;
t104 = qJDD(1) * t71;
t96 = t72 * t104;
t22 = -qJD(5) * t33 - t43 * t76 + t78 * t96;
t32 = -t45 * t76 + t78 * t99;
t23 = qJD(5) * t32 + t43 * t78 + t76 * t96;
t41 = qJD(5) + t44;
t24 = -mrSges(6,2) * t41 + mrSges(6,3) * t32;
t25 = mrSges(6,1) * t41 - mrSges(6,3) * t33;
t31 = pkin(4) * t44 - pkin(6) * t45;
t55 = t90 * qJD(1);
t95 = -g(2) * t77 + t79 * g(3);
t54 = -pkin(1) * t80 + qJDD(1) * qJ(2) + t95;
t97 = 0.2e1 * qJD(1) * qJD(2);
t93 = -t72 * g(1) + (t54 + t97) * t75;
t29 = t105 * t55 + t93;
t100 = t117 * t71 + t74 * t29;
t103 = qJDD(1) * t75;
t113 = t69 * t80;
t48 = (pkin(3) * t71 - qJ(4) * t74) * t106;
t18 = -pkin(3) * t113 - qJ(4) * t103 - t48 * t99 + t100;
t110 = t75 * t80;
t109 = -t75 * g(1) - t72 * t54;
t28 = t55 * t106 + t72 * t97 + qJDD(3) - t109;
t20 = ((-qJDD(1) * t74 - t110 * t71) * qJ(4) + (-t110 * t74 + t104) * pkin(3)) * t72 + t28;
t89 = t70 * t18 - t73 * t20;
t82 = m(6) * (-pkin(4) * t96 - pkin(6) * t102 + (t115 + t31) * t45 + t89) - t22 * mrSges(6,1) + t23 * mrSges(6,2) - t32 * t24 + t33 * t25;
t10 = m(5) * (-0.2e1 * qJD(4) * t45 - t89) - t43 * mrSges(5,3) - t45 * t30 + (mrSges(5,1) * qJDD(1) + qJD(1) * t35) * t112 - t82;
t91 = mrSges(4,1) * t71 + mrSges(4,2) * t74;
t49 = t91 * t106;
t53 = (-mrSges(4,1) * t75 - mrSges(4,3) * t111) * qJD(1);
t87 = mrSges(4,2) * t75 - mrSges(4,3) * t112;
t101 = -t44 * t115 + t18 * t73 + t70 * t20;
t14 = -pkin(4) * t102 + pkin(6) * t96 - t31 * t44 + t101;
t26 = t71 * t29;
t17 = pkin(3) * t103 - qJ(4) * t113 - t117 * t74 + t48 * t98 + qJDD(4) + t26;
t15 = (t44 * t99 - t43) * pkin(6) + (t45 * t99 + t42) * pkin(4) + t17;
t21 = -mrSges(6,1) * t32 + mrSges(6,2) * t33;
t40 = qJDD(5) + t42;
t11 = m(6) * (-t14 * t76 + t15 * t78) - t23 * mrSges(6,3) + t40 * mrSges(6,1) - t33 * t21 + t41 * t24;
t12 = m(6) * (t14 * t78 + t15 * t76) + t22 * mrSges(6,3) - t40 * mrSges(6,2) + t32 * t21 - t41 * t25;
t36 = mrSges(5,1) * t99 - mrSges(5,3) * t45;
t9 = m(5) * t101 - t42 * mrSges(5,3) - t44 * t30 + t78 * t12 - t76 * t11 + (-mrSges(5,2) * qJDD(1) - qJD(1) * t36) * t112;
t7 = m(4) * t100 + t73 * t9 - t70 * t10 + t87 * qJDD(1) + (-t112 * t49 + t53 * t75) * qJD(1);
t52 = t87 * qJD(1);
t81 = m(5) * t17 + t42 * mrSges(5,1) + t43 * mrSges(5,2) + t78 * t11 + t76 * t12 + t44 * t35 + t45 * t36;
t8 = -m(4) * t26 + (-mrSges(4,1) * qJDD(1) - qJD(1) * t52) * t75 + (m(4) * t37 + (-qJDD(1) * mrSges(4,3) + (m(4) * t118 - t49) * qJD(1)) * t72) * t74 - t81;
t4 = m(3) * t93 + t74 * t7 - t71 * t8 + (qJDD(1) * mrSges(3,3) + qJD(1) * t56) * t75;
t83 = m(4) * t28 + t73 * t10 + t70 * t9;
t88 = t52 * t71 + t53 * t74;
t6 = m(3) * t109 + ((-mrSges(3,3) - t91) * qJDD(1) + (-0.2e1 * m(3) * qJD(2) - t56 - t88) * qJD(1)) * t72 - t83;
t114 = t4 * t72 + t6 * t75;
t84 = m(3) * (-qJDD(1) * pkin(1) + t85) + t71 * t7 + t74 * t8;
t2 = m(2) * t108 + (-mrSges(2,2) + t116) * t80 + (mrSges(2,1) - t92) * qJDD(1) - t84;
t1 = m(2) * t95 - mrSges(2,1) * t80 - qJDD(1) * mrSges(2,2) + t4 * t75 - t6 * t72;
t3 = [(-m(1) - m(2)) * g(1) + t114, t1, t4, t7, t9, t12; -m(1) * g(2) + t1 * t77 + t2 * t79, t2, t6, t8, t10, t11; -m(1) * g(3) - t1 * t79 + t2 * t77, -m(2) * g(1) + t114, qJDD(1) * t92 - t116 * t80 + t84, (qJD(1) * t88 + qJDD(1) * t91) * t72 + t83, t81, t82;];
f_new = t3;
