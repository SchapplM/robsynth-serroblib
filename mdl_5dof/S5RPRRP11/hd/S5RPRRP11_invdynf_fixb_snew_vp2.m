% Calculate vector of cutting forces with Newton-Euler
% S5RPRRP11
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRRP11_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP11_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP11_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP11_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP11_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP11_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP11_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP11_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP11_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:53:29
% EndTime: 2019-12-31 18:53:32
% DurationCPUTime: 1.00s
% Computational Cost: add. (9395->150), mult. (22046->186), div. (0->0), fcn. (15364->8), ass. (0->77)
t77 = qJD(1) ^ 2;
t70 = cos(pkin(8));
t68 = t70 ^ 2;
t69 = sin(pkin(8));
t99 = t69 ^ 2 + t68;
t111 = t99 * mrSges(3,3);
t72 = sin(qJ(3));
t74 = cos(qJ(3));
t84 = t69 * t72 - t70 * t74;
t58 = t84 * qJD(1);
t85 = t69 * t74 + t70 * t72;
t59 = t85 * qJD(1);
t96 = t59 * qJD(3);
t48 = -t84 * qJDD(1) - t96;
t46 = pkin(3) * t58 - pkin(7) * t59;
t76 = qJD(3) ^ 2;
t107 = pkin(2) * t77;
t73 = sin(qJ(1));
t75 = cos(qJ(1));
t88 = -g(1) * t75 - g(2) * t73;
t60 = -pkin(1) * t77 + qJDD(1) * qJ(2) + t88;
t95 = qJD(1) * qJD(2);
t91 = -g(3) * t70 - 0.2e1 * t69 * t95;
t98 = pkin(6) * qJDD(1);
t38 = (t70 * t107 - t60 - t98) * t69 + t91;
t89 = -g(3) * t69 + (t60 + 0.2e1 * t95) * t70;
t39 = -t68 * t107 + t70 * t98 + t89;
t90 = t38 * t74 - t72 * t39;
t18 = -qJDD(3) * pkin(3) - pkin(7) * t76 + t59 * t46 - t90;
t106 = cos(qJ(4));
t97 = t58 * qJD(3);
t49 = t85 * qJDD(1) - t97;
t71 = sin(qJ(4));
t51 = t71 * qJD(3) + t106 * t59;
t24 = t51 * qJD(4) - t106 * qJDD(3) + t49 * t71;
t50 = -t106 * qJD(3) + t59 * t71;
t25 = -t50 * qJD(4) + t71 * qJDD(3) + t106 * t49;
t56 = qJD(4) + t58;
t33 = -mrSges(5,2) * t56 - mrSges(5,3) * t50;
t34 = mrSges(5,1) * t56 - mrSges(5,3) * t51;
t35 = -mrSges(6,1) * t56 + mrSges(6,2) * t51;
t32 = -mrSges(6,2) * t50 + mrSges(6,3) * t56;
t93 = m(6) * (-0.2e1 * qJD(5) * t51 + (t50 * t56 - t25) * qJ(5) + (t51 * t56 + t24) * pkin(4) + t18) + t24 * mrSges(6,1) + t50 * t32;
t110 = m(5) * t18 + t24 * mrSges(5,1) + (t34 - t35) * t51 + (mrSges(5,2) - mrSges(6,3)) * t25 + t50 * t33 + t93;
t28 = mrSges(6,1) * t50 - mrSges(6,3) * t51;
t102 = -mrSges(5,1) * t50 - mrSges(5,2) * t51 - t28;
t100 = t72 * t38 + t74 * t39;
t19 = -pkin(3) * t76 + qJDD(3) * pkin(7) - t46 * t58 + t100;
t92 = g(1) * t73 - t75 * g(2);
t87 = qJDD(2) - t92;
t78 = (-pkin(2) * t70 - pkin(1)) * qJDD(1) + (-t99 * pkin(6) - qJ(2)) * t77 + t87;
t21 = (-t49 + t97) * pkin(7) + (-t48 + t96) * pkin(3) + t78;
t103 = t106 * t19 + t71 * t21;
t104 = -mrSges(5,3) - mrSges(6,2);
t45 = qJDD(4) - t48;
t27 = pkin(4) * t50 - qJ(5) * t51;
t55 = t56 ^ 2;
t94 = m(6) * (-pkin(4) * t55 + qJ(5) * t45 + 0.2e1 * qJD(5) * t56 - t27 * t50 + t103) + t56 * t35 + t45 * mrSges(6,3);
t10 = m(5) * t103 - t45 * mrSges(5,2) + t102 * t50 + t104 * t24 - t56 * t34 + t94;
t82 = t106 * t21 - t71 * t19;
t108 = m(6) * (-t45 * pkin(4) - t55 * qJ(5) + t51 * t27 + qJDD(5) - t82);
t12 = m(5) * t82 - t108 + (t33 + t32) * t56 + t102 * t51 + (mrSges(5,1) + mrSges(6,1)) * t45 + t104 * t25;
t43 = mrSges(4,1) * t58 + mrSges(4,2) * t59;
t53 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t59;
t7 = m(4) * t100 - qJDD(3) * mrSges(4,2) + t48 * mrSges(4,3) - qJD(3) * t53 + t106 * t10 - t71 * t12 - t58 * t43;
t52 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t58;
t8 = m(4) * t90 + qJDD(3) * mrSges(4,1) - t49 * mrSges(4,3) + qJD(3) * t52 - t59 * t43 - t110;
t86 = -mrSges(3,1) * t70 + mrSges(3,2) * t69;
t83 = mrSges(3,3) * qJDD(1) + t77 * t86;
t4 = m(3) * t91 + t72 * t7 + t74 * t8 + (-m(3) * t60 - t83) * t69;
t5 = m(3) * t89 + t74 * t7 + t83 * t70 - t72 * t8;
t109 = t70 * t4 + t69 * t5;
t81 = -m(4) * t78 + t48 * mrSges(4,1) - t49 * mrSges(4,2) - t71 * t10 - t106 * t12 - t58 * t52 - t59 * t53;
t79 = m(3) * (-qJDD(1) * pkin(1) - qJ(2) * t77 + t87) - t81;
t6 = m(2) * t92 + (-mrSges(2,2) + t111) * t77 + (mrSges(2,1) - t86) * qJDD(1) - t79;
t1 = m(2) * t88 - t77 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t69 * t4 + t70 * t5;
t2 = [-m(1) * g(1) + t1 * t75 - t6 * t73, t1, t5, t7, t10, -t24 * mrSges(6,2) - t50 * t28 + t94; -m(1) * g(2) + t1 * t73 + t6 * t75, t6, t4, t8, t12, -t25 * mrSges(6,3) - t51 * t35 + t93; (-m(1) - m(2)) * g(3) + t109, -m(2) * g(3) + t109, t86 * qJDD(1) - t77 * t111 + t79, -t81, t110, -t45 * mrSges(6,1) + t25 * mrSges(6,2) + t51 * t28 - t56 * t32 + t108;];
f_new = t2;
