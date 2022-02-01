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
% m [6x1]
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
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 08:59:07
% EndTime: 2022-01-23 08:59:10
% DurationCPUTime: 1.42s
% Computational Cost: add. (10432->149), mult. (29290->220), div. (0->0), fcn. (19772->10), ass. (0->93)
t116 = -2 * qJD(3);
t70 = sin(pkin(7));
t105 = qJD(1) * t70;
t78 = qJD(1) ^ 2;
t75 = sin(qJ(1));
t77 = cos(qJ(1));
t94 = t75 * g(1) - t77 * g(2);
t81 = -t78 * qJ(2) + qJDD(2) - t94;
t73 = cos(pkin(7));
t88 = -pkin(2) * t73 - qJ(3) * t70;
t37 = (-pkin(1) + t88) * qJDD(1) + t81;
t115 = t105 * t116 + t37;
t66 = t70 ^ 2;
t67 = t73 ^ 2;
t114 = (t66 + t67) * mrSges(3,3);
t72 = cos(pkin(8));
t109 = t70 * t72;
t68 = sin(pkin(9));
t71 = cos(pkin(9));
t84 = t68 * t109 + t71 * t73;
t44 = t84 * qJD(1);
t42 = t84 * qJDD(1);
t113 = 2 * qJD(4);
t90 = -t73 * mrSges(3,1) + t70 * mrSges(3,2);
t56 = t90 * qJD(1);
t69 = sin(pkin(8));
t110 = t69 * t70;
t104 = qJD(1) * t73;
t97 = t72 * t105;
t45 = -t68 * t104 + t71 * t97;
t30 = t44 * mrSges(5,1) + t45 * mrSges(5,2);
t98 = t69 * t105;
t35 = -mrSges(5,2) * t98 - t44 * mrSges(5,3);
t43 = (t71 * t109 - t68 * t73) * qJDD(1);
t101 = t69 ^ 2 * t66 * t78;
t74 = sin(qJ(5));
t76 = cos(qJ(5));
t33 = t76 * t45 + t74 * t98;
t103 = qJDD(1) * t69;
t95 = t70 * t103;
t22 = -t33 * qJD(5) - t74 * t43 + t76 * t95;
t32 = -t74 * t45 + t76 * t98;
t23 = t32 * qJD(5) + t76 * t43 + t74 * t95;
t41 = qJD(5) + t44;
t24 = -t41 * mrSges(6,2) + t32 * mrSges(6,3);
t25 = t41 * mrSges(6,1) - t33 * mrSges(6,3);
t31 = t44 * pkin(4) - t45 * pkin(6);
t102 = qJDD(1) * t73;
t111 = t67 * t78;
t48 = (pkin(3) * t69 - qJ(4) * t72) * t105;
t55 = t88 * qJD(1);
t91 = -t77 * g(1) - t75 * g(2);
t54 = -t78 * pkin(1) + qJDD(1) * qJ(2) + t91;
t96 = 0.2e1 * qJD(1) * qJD(2);
t92 = -t70 * g(3) + (t54 + t96) * t73;
t29 = t55 * t104 + t92;
t99 = t115 * t69 + t72 * t29;
t18 = -pkin(3) * t111 - qJ(4) * t102 - t48 * t98 + t99;
t108 = t73 * t78;
t107 = -t73 * g(3) - t70 * t54;
t28 = t55 * t105 + t70 * t96 + qJDD(3) - t107;
t20 = ((-qJDD(1) * t72 - t69 * t108) * qJ(4) + (-t72 * t108 + t103) * pkin(3)) * t70 + t28;
t87 = t68 * t18 - t71 * t20;
t80 = m(6) * (-pkin(4) * t95 - pkin(6) * t101 + (t113 + t31) * t45 + t87) - t22 * mrSges(6,1) + t23 * mrSges(6,2) - t32 * t24 + t33 * t25;
t10 = m(5) * (-0.2e1 * qJD(4) * t45 - t87) - t43 * mrSges(5,3) - t45 * t30 + (mrSges(5,1) * qJDD(1) + qJD(1) * t35) * t110 - t80;
t89 = mrSges(4,1) * t69 + mrSges(4,2) * t72;
t49 = t89 * t105;
t53 = (-t73 * mrSges(4,1) - mrSges(4,3) * t109) * qJD(1);
t85 = t73 * mrSges(4,2) - mrSges(4,3) * t110;
t100 = -t44 * t113 + t71 * t18 + t68 * t20;
t14 = -pkin(4) * t101 + pkin(6) * t95 - t44 * t31 + t100;
t26 = t69 * t29;
t17 = pkin(3) * t102 - qJ(4) * t111 - t115 * t72 + t48 * t97 + qJDD(4) + t26;
t15 = (t44 * t98 - t43) * pkin(6) + (t45 * t98 + t42) * pkin(4) + t17;
t21 = -t32 * mrSges(6,1) + t33 * mrSges(6,2);
t40 = qJDD(5) + t42;
t11 = m(6) * (-t74 * t14 + t76 * t15) - t23 * mrSges(6,3) + t40 * mrSges(6,1) - t33 * t21 + t41 * t24;
t12 = m(6) * (t76 * t14 + t74 * t15) + t22 * mrSges(6,3) - t40 * mrSges(6,2) + t32 * t21 - t41 * t25;
t36 = mrSges(5,1) * t98 - t45 * mrSges(5,3);
t9 = m(5) * t100 - t42 * mrSges(5,3) - t44 * t30 + t76 * t12 - t74 * t11 + (-mrSges(5,2) * qJDD(1) - qJD(1) * t36) * t110;
t7 = m(4) * t99 + t71 * t9 - t68 * t10 + t85 * qJDD(1) + (-t49 * t110 + t73 * t53) * qJD(1);
t52 = t85 * qJD(1);
t79 = m(5) * t17 + t42 * mrSges(5,1) + t43 * mrSges(5,2) + t76 * t11 + t74 * t12 + t44 * t35 + t45 * t36;
t8 = -m(4) * t26 + (-mrSges(4,1) * qJDD(1) - qJD(1) * t52) * t73 + (m(4) * t37 + (-qJDD(1) * mrSges(4,3) + (m(4) * t116 - t49) * qJD(1)) * t70) * t72 - t79;
t4 = m(3) * t92 + t72 * t7 - t69 * t8 + (qJDD(1) * mrSges(3,3) + qJD(1) * t56) * t73;
t82 = m(4) * t28 + t71 * t10 + t68 * t9;
t86 = t52 * t69 + t53 * t72;
t6 = m(3) * t107 + ((-mrSges(3,3) - t89) * qJDD(1) + (-0.2e1 * m(3) * qJD(2) - t56 - t86) * qJD(1)) * t70 - t82;
t112 = t70 * t4 + t73 * t6;
t83 = m(3) * (-qJDD(1) * pkin(1) + t81) + t69 * t7 + t72 * t8;
t2 = m(2) * t94 + (-mrSges(2,2) + t114) * t78 + (mrSges(2,1) - t90) * qJDD(1) - t83;
t1 = m(2) * t91 - t78 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t73 * t4 - t70 * t6;
t3 = [-m(1) * g(1) + t77 * t1 - t75 * t2, t1, t4, t7, t9, t12; -m(1) * g(2) + t75 * t1 + t77 * t2, t2, t6, t8, t10, t11; (-m(1) - m(2)) * g(3) + t112, -m(2) * g(3) + t112, t90 * qJDD(1) - t78 * t114 + t83, (t86 * qJD(1) + t89 * qJDD(1)) * t70 + t82, t79, t80;];
f_new = t3;
