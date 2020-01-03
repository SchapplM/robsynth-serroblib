% Calculate vector of cutting forces with Newton-Euler
% S5RRPRP7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 20:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRPRP7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP7_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP7_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP7_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP7_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP7_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:59:58
% EndTime: 2019-12-31 20:00:01
% DurationCPUTime: 1.19s
% Computational Cost: add. (10407->162), mult. (23487->209), div. (0->0), fcn. (15404->8), ass. (0->77)
t115 = -2 * qJD(3);
t82 = qJD(1) ^ 2;
t78 = sin(qJ(1));
t80 = cos(qJ(1));
t91 = -t80 * g(1) - t78 * g(2);
t62 = -t82 * pkin(1) + qJDD(1) * pkin(6) + t91;
t77 = sin(qJ(2));
t107 = t77 * t62;
t109 = pkin(2) * t82;
t79 = cos(qJ(2));
t98 = qJD(1) * qJD(2);
t65 = t77 * qJDD(1) + t79 * t98;
t30 = qJDD(2) * pkin(2) - t65 * qJ(3) - t107 + (qJ(3) * t98 + t77 * t109 - g(3)) * t79;
t66 = t79 * qJDD(1) - t77 * t98;
t101 = qJD(1) * t77;
t67 = qJD(2) * pkin(2) - qJ(3) * t101;
t73 = t79 ^ 2;
t93 = -t77 * g(3) + t79 * t62;
t31 = t66 * qJ(3) - qJD(2) * t67 - t73 * t109 + t93;
t74 = sin(pkin(8));
t75 = cos(pkin(8));
t60 = (t74 * t79 + t75 * t77) * qJD(1);
t114 = t60 * t115 + t75 * t30 - t74 * t31;
t59 = (t74 * t77 - t75 * t79) * qJD(1);
t68 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t101;
t100 = qJD(1) * t79;
t69 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t100;
t108 = cos(qJ(4));
t76 = sin(qJ(4));
t50 = -t108 * qJD(2) + t76 * t60;
t51 = t76 * qJD(2) + t108 * t60;
t33 = t50 * mrSges(6,1) - t51 * mrSges(6,3);
t103 = -t50 * mrSges(5,1) - t51 * mrSges(5,2) - t33;
t44 = t59 * pkin(3) - t60 * pkin(7);
t81 = qJD(2) ^ 2;
t95 = t59 * t115 + t74 * t30 + t75 * t31;
t19 = -t81 * pkin(3) + qJDD(2) * pkin(7) - t59 * t44 + t95;
t48 = -t74 * t65 + t75 * t66;
t49 = t75 * t65 + t74 * t66;
t94 = t78 * g(1) - t80 * g(2);
t88 = -qJDD(1) * pkin(1) - t94;
t84 = -t66 * pkin(2) + qJDD(3) + t67 * t101 + (-qJ(3) * t73 - pkin(6)) * t82 + t88;
t21 = (qJD(2) * t59 - t49) * pkin(7) + (qJD(2) * t60 - t48) * pkin(3) + t84;
t104 = t108 * t19 + t76 * t21;
t105 = -mrSges(5,3) - mrSges(6,2);
t24 = t51 * qJD(4) - t108 * qJDD(2) + t76 * t49;
t58 = qJD(4) + t59;
t40 = t58 * mrSges(5,1) - t51 * mrSges(5,3);
t47 = qJDD(4) - t48;
t32 = t50 * pkin(4) - t51 * qJ(5);
t41 = -t58 * mrSges(6,1) + t51 * mrSges(6,2);
t57 = t58 ^ 2;
t97 = m(6) * (-t57 * pkin(4) + t47 * qJ(5) + 0.2e1 * qJD(5) * t58 - t50 * t32 + t104) + t58 * t41 + t47 * mrSges(6,3);
t10 = m(5) * t104 - t47 * mrSges(5,2) + t103 * t50 + t105 * t24 - t58 * t40 + t97;
t87 = t108 * t21 - t76 * t19;
t110 = m(6) * (-t47 * pkin(4) - t57 * qJ(5) + t51 * t32 + qJDD(5) - t87);
t25 = -t50 * qJD(4) + t76 * qJDD(2) + t108 * t49;
t38 = -t50 * mrSges(6,2) + t58 * mrSges(6,3);
t39 = -t58 * mrSges(5,2) - t50 * mrSges(5,3);
t12 = m(5) * t87 - t110 + (t39 + t38) * t58 + t103 * t51 + (mrSges(5,1) + mrSges(6,1)) * t47 + t105 * t25;
t52 = -qJD(2) * mrSges(4,2) - t59 * mrSges(4,3);
t53 = qJD(2) * mrSges(4,1) - t60 * mrSges(4,3);
t86 = -m(4) * t84 + t48 * mrSges(4,1) - t49 * mrSges(4,2) - t76 * t10 - t108 * t12 - t59 * t52 - t60 * t53;
t113 = (t77 * t68 - t79 * t69) * qJD(1) + m(3) * (-t82 * pkin(6) + t88) - t66 * mrSges(3,1) + t65 * mrSges(3,2) - t86;
t18 = -qJDD(2) * pkin(3) - t81 * pkin(7) + t60 * t44 - t114;
t96 = m(6) * (-0.2e1 * qJD(5) * t51 + (t50 * t58 - t25) * qJ(5) + (t51 * t58 + t24) * pkin(4) + t18) + t50 * t38 + t24 * mrSges(6,1);
t112 = m(5) * t18 + t24 * mrSges(5,1) + (t40 - t41) * t51 + (mrSges(5,2) - mrSges(6,3)) * t25 + t50 * t39 + t96;
t64 = (-mrSges(3,1) * t79 + mrSges(3,2) * t77) * qJD(1);
t43 = t59 * mrSges(4,1) + t60 * mrSges(4,2);
t7 = m(4) * t95 - qJDD(2) * mrSges(4,2) + t48 * mrSges(4,3) - qJD(2) * t53 + t108 * t10 - t76 * t12 - t59 * t43;
t8 = m(4) * t114 + qJDD(2) * mrSges(4,1) - t49 * mrSges(4,3) + qJD(2) * t52 - t60 * t43 - t112;
t4 = m(3) * (-t79 * g(3) - t107) - t65 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t64 * t101 + qJD(2) * t69 + t74 * t7 + t75 * t8;
t5 = m(3) * t93 - qJDD(2) * mrSges(3,2) + t66 * mrSges(3,3) - qJD(2) * t68 + t64 * t100 + t75 * t7 - t74 * t8;
t111 = t79 * t4 + t77 * t5;
t6 = m(2) * t94 + qJDD(1) * mrSges(2,1) - t82 * mrSges(2,2) - t113;
t1 = m(2) * t91 - t82 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t77 * t4 + t79 * t5;
t2 = [-m(1) * g(1) + t80 * t1 - t78 * t6, t1, t5, t7, t10, -t24 * mrSges(6,2) - t50 * t33 + t97; -m(1) * g(2) + t78 * t1 + t80 * t6, t6, t4, t8, t12, -t25 * mrSges(6,3) - t51 * t41 + t96; (-m(1) - m(2)) * g(3) + t111, -m(2) * g(3) + t111, t113, -t86, t112, -t47 * mrSges(6,1) + t25 * mrSges(6,2) + t51 * t33 - t58 * t38 + t110;];
f_new = t2;
