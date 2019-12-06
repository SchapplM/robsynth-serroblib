% Calculate vector of cutting forces with Newton-Euler
% S5RPRRP3
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
% Datum: 2019-12-05 18:04
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRRP3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP3_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP3_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP3_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP3_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:03:13
% EndTime: 2019-12-05 18:03:16
% DurationCPUTime: 0.79s
% Computational Cost: add. (6692->134), mult. (13296->172), div. (0->0), fcn. (7843->8), ass. (0->65)
t70 = sin(qJ(3));
t73 = cos(qJ(3));
t90 = qJD(1) * qJD(3);
t84 = t73 * t90;
t51 = t70 * qJDD(1) + t84;
t52 = t73 * qJDD(1) - t70 * t90;
t92 = qJD(1) * t70;
t53 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t92;
t91 = qJD(1) * t73;
t54 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t91;
t75 = qJD(1) ^ 2;
t69 = sin(qJ(4));
t72 = cos(qJ(4));
t46 = (t69 * t73 + t70 * t72) * qJD(1);
t26 = -t46 * qJD(4) - t69 * t51 + t72 * t52;
t45 = (-t69 * t70 + t72 * t73) * qJD(1);
t27 = t45 * qJD(4) + t72 * t51 + t69 * t52;
t64 = qJD(3) + qJD(4);
t36 = -t64 * mrSges(6,2) + t45 * mrSges(6,3);
t37 = -t64 * mrSges(5,2) + t45 * mrSges(5,3);
t40 = t64 * mrSges(5,1) - t46 * mrSges(5,3);
t55 = qJD(3) * pkin(3) - pkin(7) * t92;
t65 = t73 ^ 2;
t71 = sin(qJ(1));
t74 = cos(qJ(1));
t93 = t74 * g(2) + t71 * g(3);
t48 = qJDD(1) * pkin(1) + t93;
t85 = t71 * g(2) - t74 * g(3);
t50 = -t75 * pkin(1) + t85;
t67 = sin(pkin(8));
t68 = cos(pkin(8));
t81 = t68 * t48 - t67 * t50;
t78 = -qJDD(1) * pkin(2) - t81;
t76 = -t52 * pkin(3) + t55 * t92 + (-pkin(7) * t65 - pkin(6)) * t75 + t78;
t38 = t64 * pkin(4) - t46 * qJ(5);
t39 = t64 * mrSges(6,1) - t46 * mrSges(6,3);
t41 = t45 ^ 2;
t86 = m(6) * (-t26 * pkin(4) - t41 * qJ(5) + t46 * t38 + qJDD(5) + t76) + t27 * mrSges(6,2) + t46 * t39;
t77 = m(5) * t76 + t27 * mrSges(5,2) - (mrSges(5,1) + mrSges(6,1)) * t26 + t46 * t40 - (t37 + t36) * t45 + t86;
t103 = (t70 * t53 - t73 * t54) * qJD(1) + t51 * mrSges(4,2) + m(4) * (-t75 * pkin(6) + t78) - t52 * mrSges(4,1) + t77;
t94 = t67 * t48 + t68 * t50;
t32 = -t75 * pkin(2) + qJDD(1) * pkin(6) + t94;
t66 = -g(1) + qJDD(2);
t82 = -t70 * t32 + t73 * t66;
t18 = (-t51 + t84) * pkin(7) + (t70 * t73 * t75 + qJDD(3)) * pkin(3) + t82;
t96 = t73 * t32 + t70 * t66;
t19 = -t65 * t75 * pkin(3) + t52 * pkin(7) - qJD(3) * t55 + t96;
t97 = t69 * t18 + t72 * t19;
t34 = -t45 * mrSges(5,1) + t46 * mrSges(5,2);
t63 = qJDD(3) + qJDD(4);
t33 = -t45 * mrSges(6,1) + t46 * mrSges(6,2);
t87 = m(6) * (-t41 * pkin(4) + t26 * qJ(5) + 0.2e1 * qJD(5) * t45 - t64 * t38 + t97) + t45 * t33 + t26 * mrSges(6,3);
t10 = m(5) * t97 + t26 * mrSges(5,3) + t45 * t34 + (-t40 - t39) * t64 + (-mrSges(5,2) - mrSges(6,2)) * t63 + t87;
t49 = (-mrSges(4,1) * t73 + mrSges(4,2) * t70) * qJD(1);
t83 = t72 * t18 - t69 * t19;
t88 = m(6) * (-0.2e1 * qJD(5) * t46 + (t45 * t64 - t27) * qJ(5) + (t45 * t46 + t63) * pkin(4) + t83) + t64 * t36 + t63 * mrSges(6,1);
t9 = m(5) * t83 + t63 * mrSges(5,1) + t64 * t37 + (-t34 - t33) * t46 + (-mrSges(5,3) - mrSges(6,3)) * t27 + t88;
t6 = m(4) * t82 + qJDD(3) * mrSges(4,1) - t51 * mrSges(4,3) + qJD(3) * t54 + t69 * t10 - t49 * t92 + t72 * t9;
t7 = m(4) * t96 - qJDD(3) * mrSges(4,2) + t52 * mrSges(4,3) - qJD(3) * t53 + t72 * t10 + t49 * t91 - t69 * t9;
t89 = m(3) * t66 + t73 * t6 + t70 * t7;
t8 = m(3) * t81 + qJDD(1) * mrSges(3,1) - t75 * mrSges(3,2) - t103;
t3 = m(3) * t94 - t75 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t70 * t6 + t73 * t7;
t2 = m(2) * t93 + qJDD(1) * mrSges(2,1) - t75 * mrSges(2,2) + t67 * t3 + t68 * t8;
t1 = m(2) * t85 - t75 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t68 * t3 - t67 * t8;
t4 = [(-m(1) - m(2)) * g(1) + t89, t1, t3, t7, t10, -t63 * mrSges(6,2) - t64 * t39 + t87; -m(1) * g(2) - t71 * t1 - t74 * t2, t2, t8, t6, t9, -t27 * mrSges(6,3) - t46 * t33 + t88; -m(1) * g(3) + t74 * t1 - t71 * t2, -m(2) * g(1) + t89, t89, t103, t77, -t26 * mrSges(6,1) - t45 * t36 + t86;];
f_new = t4;
