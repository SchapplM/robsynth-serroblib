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
% Datum: 2019-12-05 18:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:33:42
% EndTime: 2019-12-05 18:33:45
% DurationCPUTime: 1.28s
% Computational Cost: add. (20601->129), mult. (28559->167), div. (0->0), fcn. (19171->10), ass. (0->76)
t64 = qJD(1) + qJD(2);
t58 = t64 ^ 2;
t66 = cos(pkin(9));
t62 = t66 ^ 2;
t65 = sin(pkin(9));
t92 = t65 ^ 2 + t62;
t100 = t92 * mrSges(4,3);
t99 = -m(2) - m(3);
t60 = qJDD(1) + qJDD(2);
t70 = sin(qJ(1));
t74 = cos(qJ(1));
t93 = t74 * g(2) + t70 * g(3);
t51 = qJDD(1) * pkin(1) + t93;
t75 = qJD(1) ^ 2;
t88 = t70 * g(2) - t74 * g(3);
t52 = -t75 * pkin(1) + t88;
t69 = sin(qJ(2));
t73 = cos(qJ(2));
t94 = t69 * t51 + t73 * t52;
t40 = -t58 * pkin(2) + t60 * qJ(3) + t94;
t68 = sin(qJ(4));
t72 = cos(qJ(4));
t81 = t65 * t72 + t66 * t68;
t80 = -t65 * t68 + t66 * t72;
t45 = t80 * t64;
t90 = t45 * qJD(4);
t39 = t81 * t60 + t90;
t46 = t81 * t64;
t91 = qJD(3) * t64;
t89 = -t66 * g(1) - 0.2e1 * t65 * t91;
t96 = pkin(7) * t60;
t97 = pkin(3) * t66;
t25 = (t58 * t97 - t40 - t96) * t65 + t89;
t85 = -t65 * g(1) + (t40 + 0.2e1 * t91) * t66;
t26 = -t62 * t58 * pkin(3) + t66 * t96 + t85;
t87 = t72 * t25 - t68 * t26;
t13 = (-t39 + t90) * pkin(8) + (t45 * t46 + qJDD(4)) * pkin(4) + t87;
t38 = -t46 * qJD(4) + t80 * t60;
t43 = qJD(4) * pkin(4) - t46 * pkin(8);
t44 = t45 ^ 2;
t95 = t68 * t25 + t72 * t26;
t14 = -t44 * pkin(4) + t38 * pkin(8) - qJD(4) * t43 + t95;
t67 = sin(qJ(5));
t71 = cos(qJ(5));
t31 = t71 * t45 - t67 * t46;
t19 = t31 * qJD(5) + t67 * t38 + t71 * t39;
t32 = t67 * t45 + t71 * t46;
t21 = -t31 * mrSges(6,1) + t32 * mrSges(6,2);
t63 = qJD(4) + qJD(5);
t29 = -t63 * mrSges(6,2) + t31 * mrSges(6,3);
t59 = qJDD(4) + qJDD(5);
t11 = m(6) * (t71 * t13 - t67 * t14) - t19 * mrSges(6,3) + t59 * mrSges(6,1) - t32 * t21 + t63 * t29;
t18 = -t32 * qJD(5) + t71 * t38 - t67 * t39;
t30 = t63 * mrSges(6,1) - t32 * mrSges(6,3);
t12 = m(6) * (t67 * t13 + t71 * t14) + t18 * mrSges(6,3) - t59 * mrSges(6,2) + t31 * t21 - t63 * t30;
t36 = -t45 * mrSges(5,1) + t46 * mrSges(5,2);
t41 = -qJD(4) * mrSges(5,2) + t45 * mrSges(5,3);
t8 = m(5) * t87 + qJDD(4) * mrSges(5,1) - t39 * mrSges(5,3) + qJD(4) * t41 + t71 * t11 + t67 * t12 - t46 * t36;
t83 = -t66 * mrSges(4,1) + t65 * mrSges(4,2);
t82 = t60 * mrSges(4,3) + t58 * t83;
t42 = qJD(4) * mrSges(5,1) - t46 * mrSges(5,3);
t9 = m(5) * t95 - qJDD(4) * mrSges(5,2) + t38 * mrSges(5,3) - qJD(4) * t42 - t67 * t11 + t71 * t12 + t45 * t36;
t6 = m(4) * t89 + t68 * t9 + t72 * t8 + (-m(4) * t40 - t82) * t65;
t7 = m(4) * t85 + t82 * t66 - t68 * t8 + t72 * t9;
t98 = t66 * t6 + t65 * t7;
t86 = t73 * t51 - t69 * t52;
t84 = qJDD(3) - t86;
t77 = (-pkin(2) - t97) * t60 + (-t92 * pkin(7) - qJ(3)) * t58 + t84;
t79 = t18 * mrSges(6,1) + t31 * t29 - m(6) * (-t38 * pkin(4) - t44 * pkin(8) + t46 * t43 + t77) - t19 * mrSges(6,2) - t32 * t30;
t78 = -m(5) * t77 + t38 * mrSges(5,1) - t39 * mrSges(5,2) + t45 * t41 - t46 * t42 + t79;
t76 = m(4) * (-t60 * pkin(2) - t58 * qJ(3) + t84) - t78;
t10 = -t76 + (-mrSges(3,2) + t100) * t58 + (mrSges(3,1) - t83) * t60 + m(3) * t86;
t3 = m(3) * t94 - t58 * mrSges(3,1) - t60 * mrSges(3,2) - t65 * t6 + t66 * t7;
t2 = m(2) * t93 + qJDD(1) * mrSges(2,1) - t75 * mrSges(2,2) + t73 * t10 + t69 * t3;
t1 = m(2) * t88 - t75 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t69 * t10 + t73 * t3;
t4 = [(-m(1) + t99) * g(1) + t98, t1, t3, t7, t9, t12; -m(1) * g(2) - t70 * t1 - t74 * t2, t2, t10, t6, t8, t11; -m(1) * g(3) + t74 * t1 - t70 * t2, t99 * g(1) + t98, -m(3) * g(1) + t98, -t58 * t100 + t83 * t60 + t76, -t78, -t79;];
f_new = t4;
