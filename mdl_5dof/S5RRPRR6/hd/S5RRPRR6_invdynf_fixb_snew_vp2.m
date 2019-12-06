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
% Datum: 2019-12-05 18:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:35:46
% EndTime: 2019-12-05 18:35:48
% DurationCPUTime: 1.01s
% Computational Cost: add. (14930->128), mult. (20438->176), div. (0->0), fcn. (12529->10), ass. (0->80)
t66 = sin(pkin(9));
t62 = t66 ^ 2;
t67 = cos(pkin(9));
t103 = (t67 ^ 2 + t62) * mrSges(4,3);
t102 = -m(2) - m(3);
t61 = qJDD(1) + qJDD(2);
t100 = t61 * mrSges(4,3);
t64 = qJD(1) + qJD(2);
t82 = -t67 * mrSges(4,1) + t66 * mrSges(4,2);
t43 = t82 * t64;
t83 = -pkin(3) * t67 - pkin(7) * t66;
t44 = t83 * t64;
t87 = 0.2e1 * qJD(3) * t64;
t60 = t64 ^ 2;
t71 = sin(qJ(1));
t75 = cos(qJ(1));
t93 = t75 * g(2) + t71 * g(3);
t48 = qJDD(1) * pkin(1) + t93;
t76 = qJD(1) ^ 2;
t86 = t71 * g(2) - t75 * g(3);
t49 = -t76 * pkin(1) + t86;
t70 = sin(qJ(2));
t74 = cos(qJ(2));
t94 = t70 * t48 + t74 * t49;
t31 = -t60 * pkin(2) + t61 * qJ(3) + t94;
t95 = -t67 * g(1) - t66 * t31;
t98 = t64 * t66;
t18 = t44 * t98 + t66 * t87 - t95;
t69 = sin(qJ(4));
t73 = cos(qJ(4));
t91 = qJD(4) * t64;
t41 = (-t61 * t69 - t73 * t91) * t66;
t42 = (t61 * t73 - t69 * t91) * t66;
t68 = sin(qJ(5));
t72 = cos(qJ(5));
t35 = (-t68 * t69 + t72 * t73) * t98;
t21 = -t35 * qJD(5) + t72 * t41 - t68 * t42;
t34 = (-t68 * t73 - t69 * t72) * t98;
t22 = t34 * qJD(5) + t68 * t41 + t72 * t42;
t97 = t67 * t64;
t53 = qJD(4) - t97;
t51 = qJD(5) + t53;
t32 = -t51 * mrSges(6,2) + t34 * mrSges(6,3);
t33 = t51 * mrSges(6,1) - t35 * mrSges(6,3);
t88 = t73 * t98;
t40 = t53 * pkin(4) - pkin(8) * t88;
t99 = t62 * t60;
t90 = t69 ^ 2 * t99;
t78 = t21 * mrSges(6,1) + t34 * t32 - m(6) * (-t41 * pkin(4) - pkin(8) * t90 + t40 * t88 + t18) - t22 * mrSges(6,2) - t35 * t33;
t77 = m(5) * t18 - t41 * mrSges(5,1) + t42 * mrSges(5,2) - t78;
t89 = t69 * t98;
t37 = -t53 * mrSges(5,2) - mrSges(5,3) * t89;
t38 = t53 * mrSges(5,1) - mrSges(5,3) * t88;
t81 = t37 * t69 + t38 * t73;
t10 = m(4) * t95 + (-t100 + (-0.2e1 * m(4) * qJD(3) - t43 - t81) * t64) * t66 - t77;
t84 = -t66 * g(1) + (t31 + t87) * t67;
t19 = t44 * t97 + t84;
t85 = t74 * t48 - t70 * t49;
t79 = -t60 * qJ(3) + qJDD(3) - t85;
t26 = (-pkin(2) + t83) * t61 + t79;
t25 = t73 * t26;
t52 = -t67 * t61 + qJDD(4);
t13 = t52 * pkin(4) - t42 * pkin(8) + t25 + (-pkin(4) * t73 * t99 - pkin(8) * t53 * t98 - t19) * t69;
t96 = t73 * t19 + t69 * t26;
t14 = -pkin(4) * t90 + t41 * pkin(8) - t53 * t40 + t96;
t27 = -t34 * mrSges(6,1) + t35 * mrSges(6,2);
t50 = qJDD(5) + t52;
t11 = m(6) * (t72 * t13 - t68 * t14) - t22 * mrSges(6,3) + t50 * mrSges(6,1) - t35 * t27 + t51 * t32;
t12 = m(6) * (t68 * t13 + t72 * t14) + t21 * mrSges(6,3) - t50 * mrSges(6,2) + t34 * t27 - t51 * t33;
t39 = (mrSges(5,1) * t69 + mrSges(5,2) * t73) * t98;
t7 = m(5) * (-t69 * t19 + t25) - t42 * mrSges(5,3) + t52 * mrSges(5,1) - t39 * t88 + t53 * t37 + t68 * t12 + t72 * t11;
t8 = m(5) * t96 - t52 * mrSges(5,2) + t41 * mrSges(5,3) - t68 * t11 + t72 * t12 - t53 * t38 - t39 * t89;
t6 = m(4) * t84 + t73 * t8 - t69 * t7 + (t64 * t43 + t100) * t67;
t101 = t67 * t10 + t66 * t6;
t80 = m(4) * (-t61 * pkin(2) + t79) + t69 * t8 + t73 * t7;
t4 = m(3) * t85 + (mrSges(3,1) - t82) * t61 + (-mrSges(3,2) + t103) * t60 - t80;
t3 = m(3) * t94 - t60 * mrSges(3,1) - t61 * mrSges(3,2) - t66 * t10 + t67 * t6;
t2 = m(2) * t93 + qJDD(1) * mrSges(2,1) - t76 * mrSges(2,2) + t70 * t3 + t74 * t4;
t1 = m(2) * t86 - t76 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t74 * t3 - t70 * t4;
t5 = [(-m(1) + t102) * g(1) + t101, t1, t3, t6, t8, t12; -m(1) * g(2) - t71 * t1 - t75 * t2, t2, t4, t10, t7, t11; -m(1) * g(3) + t75 * t1 - t71 * t2, t102 * g(1) + t101, -m(3) * g(1) + t101, -t60 * t103 + t82 * t61 + t80, t81 * t98 + t77, -t78;];
f_new = t5;
