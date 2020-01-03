% Calculate vector of cutting forces with Newton-Euler
% S5RRPRR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2020-01-03 12:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRPRR5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR5_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR5_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR5_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR5_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:03:36
% EndTime: 2020-01-03 12:03:39
% DurationCPUTime: 1.27s
% Computational Cost: add. (20601->129), mult. (28559->167), div. (0->0), fcn. (19171->10), ass. (0->76)
t62 = qJD(1) + qJD(2);
t56 = t62 ^ 2;
t64 = cos(pkin(9));
t60 = t64 ^ 2;
t63 = sin(pkin(9));
t91 = t63 ^ 2 + t60;
t98 = t91 * mrSges(4,3);
t97 = -m(2) - m(3);
t58 = qJDD(1) + qJDD(2);
t68 = sin(qJ(1));
t72 = cos(qJ(1));
t83 = -t72 * g(2) - t68 * g(3);
t51 = qJDD(1) * pkin(1) + t83;
t73 = qJD(1) ^ 2;
t87 = -t68 * g(2) + t72 * g(3);
t52 = -t73 * pkin(1) + t87;
t67 = sin(qJ(2));
t71 = cos(qJ(2));
t92 = t67 * t51 + t71 * t52;
t40 = -t56 * pkin(2) + t58 * qJ(3) + t92;
t66 = sin(qJ(4));
t70 = cos(qJ(4));
t79 = t63 * t70 + t64 * t66;
t78 = -t63 * t66 + t64 * t70;
t45 = t78 * t62;
t89 = t45 * qJD(4);
t39 = t79 * t58 + t89;
t46 = t79 * t62;
t90 = qJD(3) * t62;
t88 = -t64 * g(1) - 0.2e1 * t63 * t90;
t94 = pkin(7) * t58;
t95 = pkin(3) * t64;
t25 = (t56 * t95 - t40 - t94) * t63 + t88;
t84 = -t63 * g(1) + (t40 + 0.2e1 * t90) * t64;
t26 = -t60 * t56 * pkin(3) + t64 * t94 + t84;
t86 = t70 * t25 - t66 * t26;
t13 = (-t39 + t89) * pkin(8) + (t45 * t46 + qJDD(4)) * pkin(4) + t86;
t38 = -t46 * qJD(4) + t78 * t58;
t43 = qJD(4) * pkin(4) - t46 * pkin(8);
t44 = t45 ^ 2;
t93 = t66 * t25 + t70 * t26;
t14 = -t44 * pkin(4) + t38 * pkin(8) - qJD(4) * t43 + t93;
t65 = sin(qJ(5));
t69 = cos(qJ(5));
t31 = t69 * t45 - t65 * t46;
t19 = t31 * qJD(5) + t65 * t38 + t69 * t39;
t32 = t65 * t45 + t69 * t46;
t21 = -t31 * mrSges(6,1) + t32 * mrSges(6,2);
t61 = qJD(4) + qJD(5);
t29 = -t61 * mrSges(6,2) + t31 * mrSges(6,3);
t57 = qJDD(4) + qJDD(5);
t11 = m(6) * (t69 * t13 - t65 * t14) - t19 * mrSges(6,3) + t57 * mrSges(6,1) - t32 * t21 + t61 * t29;
t18 = -t32 * qJD(5) + t69 * t38 - t65 * t39;
t30 = t61 * mrSges(6,1) - t32 * mrSges(6,3);
t12 = m(6) * (t65 * t13 + t69 * t14) + t18 * mrSges(6,3) - t57 * mrSges(6,2) + t31 * t21 - t61 * t30;
t36 = -t45 * mrSges(5,1) + t46 * mrSges(5,2);
t41 = -qJD(4) * mrSges(5,2) + t45 * mrSges(5,3);
t8 = m(5) * t86 + qJDD(4) * mrSges(5,1) - t39 * mrSges(5,3) + qJD(4) * t41 + t69 * t11 + t65 * t12 - t46 * t36;
t81 = -t64 * mrSges(4,1) + t63 * mrSges(4,2);
t80 = t58 * mrSges(4,3) + t56 * t81;
t42 = qJD(4) * mrSges(5,1) - t46 * mrSges(5,3);
t9 = m(5) * t93 - qJDD(4) * mrSges(5,2) + t38 * mrSges(5,3) - qJD(4) * t42 - t65 * t11 + t69 * t12 + t45 * t36;
t6 = m(4) * t88 + t66 * t9 + t70 * t8 + (-m(4) * t40 - t80) * t63;
t7 = m(4) * t84 + t80 * t64 - t66 * t8 + t70 * t9;
t96 = t64 * t6 + t63 * t7;
t85 = t71 * t51 - t67 * t52;
t82 = qJDD(3) - t85;
t75 = (-pkin(2) - t95) * t58 + (-t91 * pkin(7) - qJ(3)) * t56 + t82;
t77 = t18 * mrSges(6,1) + t31 * t29 - m(6) * (-t38 * pkin(4) - t44 * pkin(8) + t46 * t43 + t75) - t19 * mrSges(6,2) - t32 * t30;
t76 = -m(5) * t75 + t38 * mrSges(5,1) - t39 * mrSges(5,2) + t45 * t41 - t46 * t42 + t77;
t74 = m(4) * (-t58 * pkin(2) - t56 * qJ(3) + t82) - t76;
t10 = -t74 + (-mrSges(3,2) + t98) * t56 + (mrSges(3,1) - t81) * t58 + m(3) * t85;
t3 = m(3) * t92 - t56 * mrSges(3,1) - t58 * mrSges(3,2) - t63 * t6 + t64 * t7;
t2 = m(2) * t83 + qJDD(1) * mrSges(2,1) - t73 * mrSges(2,2) + t71 * t10 + t67 * t3;
t1 = m(2) * t87 - t73 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t67 * t10 + t71 * t3;
t4 = [(-m(1) + t97) * g(1) + t96, t1, t3, t7, t9, t12; -m(1) * g(2) + t68 * t1 + t72 * t2, t2, t10, t6, t8, t11; -m(1) * g(3) - t72 * t1 + t68 * t2, t97 * g(1) + t96, -m(3) * g(1) + t96, -t56 * t98 + t81 * t58 + t74, -t76, -t77;];
f_new = t4;
