% Calculate vector of cutting forces with Newton-Euler
% S5RPRRP9
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
% Datum: 2019-12-31 18:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRRP9_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP9_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP9_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP9_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP9_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP9_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP9_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP9_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP9_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:48:39
% EndTime: 2019-12-31 18:48:42
% DurationCPUTime: 1.22s
% Computational Cost: add. (12129->149), mult. (29250->187), div. (0->0), fcn. (20740->8), ass. (0->76)
t78 = qJD(1) ^ 2;
t72 = cos(pkin(8));
t69 = t72 ^ 2;
t71 = sin(pkin(8));
t100 = t71 ^ 2 + t69;
t109 = t100 * mrSges(3,3);
t105 = cos(qJ(4));
t74 = sin(qJ(3));
t76 = cos(qJ(3));
t87 = -t71 * t74 + t72 * t76;
t58 = t87 * qJD(1);
t88 = t71 * t76 + t72 * t74;
t59 = t88 * qJD(1);
t73 = sin(qJ(4));
t40 = -t105 * t58 + t73 * t59;
t41 = t105 * t59 + t73 * t58;
t28 = t40 * mrSges(6,1) - t41 * mrSges(6,3);
t102 = -t40 * mrSges(5,1) - t41 * mrSges(5,2) - t28;
t104 = -mrSges(5,3) - mrSges(6,2);
t27 = t40 * pkin(4) - t41 * qJ(5);
t70 = qJD(3) + qJD(4);
t66 = t70 ^ 2;
t67 = qJDD(3) + qJDD(4);
t98 = t58 * qJD(3);
t50 = t88 * qJDD(1) + t98;
t106 = pkin(2) * t78;
t75 = sin(qJ(1));
t77 = cos(qJ(1));
t91 = -t77 * g(1) - t75 * g(2);
t60 = -t78 * pkin(1) + qJDD(1) * qJ(2) + t91;
t97 = qJD(1) * qJD(2);
t94 = -t72 * g(3) - 0.2e1 * t71 * t97;
t99 = pkin(6) * qJDD(1);
t38 = (t72 * t106 - t60 - t99) * t71 + t94;
t92 = -t71 * g(3) + (t60 + 0.2e1 * t97) * t72;
t39 = -t69 * t106 + t72 * t99 + t92;
t93 = t76 * t38 - t74 * t39;
t16 = (-t50 + t98) * pkin(7) + (t58 * t59 + qJDD(3)) * pkin(3) + t93;
t101 = t74 * t38 + t76 * t39;
t49 = -t59 * qJD(3) + t87 * qJDD(1);
t53 = qJD(3) * pkin(3) - t59 * pkin(7);
t57 = t58 ^ 2;
t18 = -t57 * pkin(3) + t49 * pkin(7) - qJD(3) * t53 + t101;
t85 = t105 * t16 - t73 * t18;
t107 = m(6) * (-t67 * pkin(4) - t66 * qJ(5) + t41 * t27 + qJDD(5) - t85);
t24 = -t40 * qJD(4) + t105 * t50 + t73 * t49;
t33 = -t70 * mrSges(5,2) - t40 * mrSges(5,3);
t36 = -t40 * mrSges(6,2) + t70 * mrSges(6,3);
t10 = m(5) * t85 - t107 + (t33 + t36) * t70 + (mrSges(5,1) + mrSges(6,1)) * t67 + t102 * t41 + t104 * t24;
t46 = -t58 * mrSges(4,1) + t59 * mrSges(4,2);
t51 = -qJD(3) * mrSges(4,2) + t58 * mrSges(4,3);
t103 = t105 * t18 + t73 * t16;
t23 = t41 * qJD(4) - t105 * t49 + t73 * t50;
t34 = t70 * mrSges(5,1) - t41 * mrSges(5,3);
t35 = -t70 * mrSges(6,1) + t41 * mrSges(6,2);
t96 = m(6) * (-t66 * pkin(4) + t67 * qJ(5) + 0.2e1 * qJD(5) * t70 - t40 * t27 + t103) + t70 * t35 + t67 * mrSges(6,3);
t9 = m(5) * t103 - t67 * mrSges(5,2) + t102 * t40 + t104 * t23 - t70 * t34 + t96;
t6 = m(4) * t93 + qJDD(3) * mrSges(4,1) - t50 * mrSges(4,3) + qJD(3) * t51 + t105 * t10 - t59 * t46 + t73 * t9;
t52 = qJD(3) * mrSges(4,1) - t59 * mrSges(4,3);
t7 = m(4) * t101 - qJDD(3) * mrSges(4,2) + t49 * mrSges(4,3) - qJD(3) * t52 - t73 * t10 + t105 * t9 + t58 * t46;
t89 = -t72 * mrSges(3,1) + t71 * mrSges(3,2);
t86 = qJDD(1) * mrSges(3,3) + t78 * t89;
t4 = m(3) * t94 + t74 * t7 + t76 * t6 + (-m(3) * t60 - t86) * t71;
t5 = m(3) * t92 - t74 * t6 + t76 * t7 + t86 * t72;
t108 = t72 * t4 + t71 * t5;
t95 = t75 * g(1) - t77 * g(2);
t90 = qJDD(2) - t95;
t82 = (-pkin(2) * t72 - pkin(1)) * qJDD(1) + (-t100 * pkin(6) - qJ(2)) * t78 + t90;
t81 = -t49 * pkin(3) - t57 * pkin(7) + t59 * t53 + t82;
t84 = t24 * mrSges(6,3) + t41 * t35 - m(6) * (-0.2e1 * qJD(5) * t41 + (t40 * t70 - t24) * qJ(5) + (t41 * t70 + t23) * pkin(4) + t81) - t23 * mrSges(6,1) - t40 * t36;
t83 = m(5) * t81 + t23 * mrSges(5,1) + t24 * mrSges(5,2) + t40 * t33 + t41 * t34 - t84;
t80 = -m(4) * t82 + t49 * mrSges(4,1) - t50 * mrSges(4,2) + t58 * t51 - t59 * t52 - t83;
t79 = m(3) * (-qJDD(1) * pkin(1) - t78 * qJ(2) + t90) - t80;
t8 = (-mrSges(2,2) + t109) * t78 + (mrSges(2,1) - t89) * qJDD(1) + m(2) * t95 - t79;
t1 = m(2) * t91 - t78 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t71 * t4 + t72 * t5;
t2 = [-m(1) * g(1) + t77 * t1 - t75 * t8, t1, t5, t7, t9, -t23 * mrSges(6,2) - t40 * t28 + t96; -m(1) * g(2) + t75 * t1 + t77 * t8, t8, t4, t6, t10, -t84; (-m(1) - m(2)) * g(3) + t108, -m(2) * g(3) + t108, t89 * qJDD(1) - t78 * t109 + t79, -t80, t83, -t67 * mrSges(6,1) + t24 * mrSges(6,2) + t41 * t28 - t70 * t36 + t107;];
f_new = t2;
