% Calculate vector of cutting forces with Newton-Euler
% S5RRRPR3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-05 18:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRPR3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR3_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR3_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR3_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR3_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:42:23
% EndTime: 2019-12-05 18:42:27
% DurationCPUTime: 1.42s
% Computational Cost: add. (23161->141), mult. (31070->190), div. (0->0), fcn. (19687->10), ass. (0->76)
t66 = qJDD(1) + qJDD(2);
t73 = sin(qJ(3));
t77 = cos(qJ(3));
t68 = qJD(1) + qJD(2);
t92 = qJD(3) * t68;
t51 = t73 * t66 + t77 * t92;
t52 = t77 * t66 - t73 * t92;
t97 = t68 * t73;
t59 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t97;
t96 = t68 * t77;
t60 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t96;
t64 = t68 ^ 2;
t70 = sin(pkin(9));
t71 = cos(pkin(9));
t39 = -t70 * t51 + t71 * t52;
t40 = t71 * t51 + t70 * t52;
t47 = (-t70 * t73 + t71 * t77) * t68;
t41 = -qJD(3) * mrSges(5,2) + t47 * mrSges(5,3);
t48 = (t70 * t77 + t71 * t73) * t68;
t42 = qJD(3) * mrSges(5,1) - t48 * mrSges(5,3);
t58 = qJD(3) * pkin(3) - qJ(4) * t97;
t69 = t77 ^ 2;
t75 = sin(qJ(1));
t79 = cos(qJ(1));
t93 = t79 * g(2) + t75 * g(3);
t56 = qJDD(1) * pkin(1) + t93;
t80 = qJD(1) ^ 2;
t89 = t75 * g(2) - t79 * g(3);
t57 = -t80 * pkin(1) + t89;
t74 = sin(qJ(2));
t78 = cos(qJ(2));
t88 = t78 * t56 - t74 * t57;
t85 = -t66 * pkin(2) - t88;
t82 = -t52 * pkin(3) + qJDD(4) + t58 * t97 + (-qJ(4) * t69 - pkin(7)) * t64 + t85;
t72 = sin(qJ(5));
t76 = cos(qJ(5));
t32 = t72 * t47 + t76 * t48;
t18 = -t32 * qJD(5) + t76 * t39 - t72 * t40;
t31 = t76 * t47 - t72 * t48;
t19 = t31 * qJD(5) + t72 * t39 + t76 * t40;
t67 = qJD(3) + qJD(5);
t29 = -t67 * mrSges(6,2) + t31 * mrSges(6,3);
t30 = t67 * mrSges(6,1) - t32 * mrSges(6,3);
t43 = qJD(3) * pkin(4) - t48 * pkin(8);
t46 = t47 ^ 2;
t84 = t18 * mrSges(6,1) + t31 * t29 - m(6) * (-t39 * pkin(4) - t46 * pkin(8) + t48 * t43 + t82) - t19 * mrSges(6,2) - t32 * t30;
t83 = -m(5) * t82 + t39 * mrSges(5,1) - t40 * mrSges(5,2) + t47 * t41 - t48 * t42 + t84;
t101 = (t73 * t59 - t77 * t60) * t68 + m(4) * (-t64 * pkin(7) + t85) - t52 * mrSges(4,1) + t51 * mrSges(4,2) - t83;
t100 = -m(2) - m(3);
t50 = (-mrSges(4,1) * t77 + mrSges(4,2) * t73) * t68;
t94 = t74 * t56 + t78 * t57;
t37 = -t64 * pkin(2) + t66 * pkin(7) + t94;
t95 = t73 * t37;
t98 = pkin(3) * t64;
t25 = qJDD(3) * pkin(3) - t51 * qJ(4) - t95 + (qJ(4) * t92 + t73 * t98 - g(1)) * t77;
t90 = -t73 * g(1) + t77 * t37;
t26 = t52 * qJ(4) - qJD(3) * t58 - t69 * t98 + t90;
t87 = -0.2e1 * qJD(4) * t48 + t71 * t25 - t70 * t26;
t13 = (qJD(3) * t47 - t40) * pkin(8) + (t47 * t48 + qJDD(3)) * pkin(4) + t87;
t91 = 0.2e1 * qJD(4) * t47 + t70 * t25 + t71 * t26;
t14 = -t46 * pkin(4) + t39 * pkin(8) - qJD(3) * t43 + t91;
t21 = -t31 * mrSges(6,1) + t32 * mrSges(6,2);
t65 = qJDD(3) + qJDD(5);
t11 = m(6) * (t76 * t13 - t72 * t14) - t19 * mrSges(6,3) + t65 * mrSges(6,1) - t32 * t21 + t67 * t29;
t12 = m(6) * (t72 * t13 + t76 * t14) + t18 * mrSges(6,3) - t65 * mrSges(6,2) + t31 * t21 - t67 * t30;
t35 = -t47 * mrSges(5,1) + t48 * mrSges(5,2);
t8 = m(5) * t87 + qJDD(3) * mrSges(5,1) - t40 * mrSges(5,3) + qJD(3) * t41 + t76 * t11 + t72 * t12 - t48 * t35;
t9 = m(5) * t91 - qJDD(3) * mrSges(5,2) + t39 * mrSges(5,3) - qJD(3) * t42 - t72 * t11 + t76 * t12 + t47 * t35;
t6 = m(4) * (-t77 * g(1) - t95) - t51 * mrSges(4,3) + qJDD(3) * mrSges(4,1) - t50 * t97 + qJD(3) * t60 + t70 * t9 + t71 * t8;
t7 = m(4) * t90 - qJDD(3) * mrSges(4,2) + t52 * mrSges(4,3) - qJD(3) * t59 + t50 * t96 - t70 * t8 + t71 * t9;
t99 = t77 * t6 + t73 * t7;
t10 = m(3) * t88 + t66 * mrSges(3,1) - t64 * mrSges(3,2) - t101;
t3 = m(3) * t94 - t64 * mrSges(3,1) - t66 * mrSges(3,2) - t73 * t6 + t77 * t7;
t2 = m(2) * t93 + qJDD(1) * mrSges(2,1) - t80 * mrSges(2,2) + t78 * t10 + t74 * t3;
t1 = m(2) * t89 - t80 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t74 * t10 + t78 * t3;
t4 = [(-m(1) + t100) * g(1) + t99, t1, t3, t7, t9, t12; -m(1) * g(2) - t75 * t1 - t79 * t2, t2, t10, t6, t8, t11; -m(1) * g(3) + t79 * t1 - t75 * t2, t100 * g(1) + t99, -m(3) * g(1) + t99, t101, -t83, -t84;];
f_new = t4;
