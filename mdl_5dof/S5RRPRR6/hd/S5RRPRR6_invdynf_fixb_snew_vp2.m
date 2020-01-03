% Calculate vector of cutting forces with Newton-Euler
% S5RRPRR6
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
% Datum: 2020-01-03 12:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRPRR6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR6_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR6_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR6_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR6_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:05:34
% EndTime: 2020-01-03 12:05:36
% DurationCPUTime: 1.01s
% Computational Cost: add. (14930->128), mult. (20438->176), div. (0->0), fcn. (12529->10), ass. (0->80)
t64 = sin(pkin(9));
t60 = t64 ^ 2;
t65 = cos(pkin(9));
t101 = (t65 ^ 2 + t60) * mrSges(4,3);
t100 = -m(2) - m(3);
t62 = qJD(1) + qJD(2);
t80 = -t65 * mrSges(4,1) + t64 * mrSges(4,2);
t43 = t80 * t62;
t82 = -pkin(3) * t65 - pkin(7) * t64;
t44 = t82 * t62;
t86 = 0.2e1 * qJD(3) * t62;
t58 = t62 ^ 2;
t59 = qJDD(1) + qJDD(2);
t69 = sin(qJ(1));
t73 = cos(qJ(1));
t81 = -t73 * g(2) - t69 * g(3);
t48 = qJDD(1) * pkin(1) + t81;
t74 = qJD(1) ^ 2;
t85 = -t69 * g(2) + t73 * g(3);
t49 = -t74 * pkin(1) + t85;
t68 = sin(qJ(2));
t72 = cos(qJ(2));
t92 = t68 * t48 + t72 * t49;
t31 = -t58 * pkin(2) + t59 * qJ(3) + t92;
t93 = -t65 * g(1) - t64 * t31;
t96 = t62 * t64;
t18 = t44 * t96 + t64 * t86 - t93;
t67 = sin(qJ(4));
t71 = cos(qJ(4));
t90 = qJD(4) * t62;
t41 = (-t59 * t67 - t71 * t90) * t64;
t42 = (t59 * t71 - t67 * t90) * t64;
t66 = sin(qJ(5));
t70 = cos(qJ(5));
t35 = (-t66 * t67 + t70 * t71) * t96;
t21 = -t35 * qJD(5) + t70 * t41 - t66 * t42;
t34 = (-t66 * t71 - t67 * t70) * t96;
t22 = t34 * qJD(5) + t66 * t41 + t70 * t42;
t95 = t65 * t62;
t53 = qJD(4) - t95;
t51 = qJD(5) + t53;
t32 = -t51 * mrSges(6,2) + t34 * mrSges(6,3);
t33 = t51 * mrSges(6,1) - t35 * mrSges(6,3);
t87 = t71 * t96;
t40 = t53 * pkin(4) - pkin(8) * t87;
t97 = t60 * t58;
t89 = t67 ^ 2 * t97;
t76 = t21 * mrSges(6,1) + t34 * t32 - m(6) * (-t41 * pkin(4) - pkin(8) * t89 + t40 * t87 + t18) - t22 * mrSges(6,2) - t35 * t33;
t75 = m(5) * t18 - t41 * mrSges(5,1) + t42 * mrSges(5,2) - t76;
t88 = t67 * t96;
t37 = -t53 * mrSges(5,2) - mrSges(5,3) * t88;
t38 = t53 * mrSges(5,1) - mrSges(5,3) * t87;
t79 = t37 * t67 + t38 * t71;
t98 = t59 * mrSges(4,3);
t10 = m(4) * t93 + (-t98 + (-0.2e1 * m(4) * qJD(3) - t43 - t79) * t62) * t64 - t75;
t83 = -t64 * g(1) + (t31 + t86) * t65;
t19 = t44 * t95 + t83;
t84 = t72 * t48 - t68 * t49;
t77 = -t58 * qJ(3) + qJDD(3) - t84;
t26 = (-pkin(2) + t82) * t59 + t77;
t25 = t71 * t26;
t52 = -t65 * t59 + qJDD(4);
t13 = t52 * pkin(4) - t42 * pkin(8) + t25 + (-pkin(4) * t71 * t97 - pkin(8) * t53 * t96 - t19) * t67;
t94 = t71 * t19 + t67 * t26;
t14 = -pkin(4) * t89 + t41 * pkin(8) - t53 * t40 + t94;
t27 = -t34 * mrSges(6,1) + t35 * mrSges(6,2);
t50 = qJDD(5) + t52;
t11 = m(6) * (t70 * t13 - t66 * t14) - t22 * mrSges(6,3) + t50 * mrSges(6,1) - t35 * t27 + t51 * t32;
t12 = m(6) * (t66 * t13 + t70 * t14) + t21 * mrSges(6,3) - t50 * mrSges(6,2) + t34 * t27 - t51 * t33;
t39 = (mrSges(5,1) * t67 + mrSges(5,2) * t71) * t96;
t7 = m(5) * (-t67 * t19 + t25) - t42 * mrSges(5,3) + t52 * mrSges(5,1) - t39 * t87 + t53 * t37 + t66 * t12 + t70 * t11;
t8 = m(5) * t94 - t52 * mrSges(5,2) + t41 * mrSges(5,3) - t66 * t11 + t70 * t12 - t53 * t38 - t39 * t88;
t6 = m(4) * t83 + t71 * t8 - t67 * t7 + (t62 * t43 + t98) * t65;
t99 = t65 * t10 + t64 * t6;
t78 = m(4) * (-t59 * pkin(2) + t77) + t67 * t8 + t71 * t7;
t4 = m(3) * t84 + (mrSges(3,1) - t80) * t59 + (-mrSges(3,2) + t101) * t58 - t78;
t3 = m(3) * t92 - t58 * mrSges(3,1) - t59 * mrSges(3,2) - t64 * t10 + t65 * t6;
t2 = m(2) * t81 + qJDD(1) * mrSges(2,1) - t74 * mrSges(2,2) + t68 * t3 + t72 * t4;
t1 = m(2) * t85 - t74 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t72 * t3 - t68 * t4;
t5 = [(-m(1) + t100) * g(1) + t99, t1, t3, t6, t8, t12; -m(1) * g(2) + t69 * t1 + t73 * t2, t2, t4, t10, t7, t11; -m(1) * g(3) - t73 * t1 + t69 * t2, t100 * g(1) + t99, -m(3) * g(1) + t99, -t58 * t101 + t80 * t59 + t78, t79 * t96 + t75, -t76;];
f_new = t5;
