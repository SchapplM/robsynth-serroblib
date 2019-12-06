% Calculate vector of cutting forces with Newton-Euler
% S5PRRRP7
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRRRP7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP7_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP7_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP7_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP7_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP7_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:54:12
% EndTime: 2019-12-05 16:54:14
% DurationCPUTime: 0.82s
% Computational Cost: add. (8551->132), mult. (16240->172), div. (0->0), fcn. (10573->10), ass. (0->71)
t64 = sin(pkin(9));
t66 = cos(pkin(9));
t57 = -t66 * g(1) - t64 * g(2);
t63 = -g(3) + qJDD(1);
t65 = sin(pkin(5));
t70 = sin(qJ(2));
t73 = cos(qJ(2));
t56 = t64 * g(1) - t66 * g(2);
t67 = cos(pkin(5));
t97 = t56 * t67;
t101 = -t70 * t57 + (t63 * t65 + t97) * t73;
t68 = sin(qJ(4));
t71 = cos(qJ(4));
t69 = sin(qJ(3));
t91 = qJD(2) * t69;
t50 = t71 * qJD(3) - t68 * t91;
t72 = cos(qJ(3));
t89 = qJD(2) * qJD(3);
t82 = t72 * t89;
t54 = t69 * qJDD(2) + t82;
t32 = t50 * qJD(4) + t68 * qJDD(3) + t71 * t54;
t51 = t68 * qJD(3) + t71 * t91;
t34 = -t50 * mrSges(6,1) + t51 * mrSges(6,2);
t35 = -t50 * mrSges(5,1) + t51 * mrSges(5,2);
t90 = t72 * qJD(2);
t61 = qJD(4) - t90;
t38 = -t61 * mrSges(5,2) + t50 * mrSges(5,3);
t83 = t69 * t89;
t55 = t72 * qJDD(2) - t83;
t47 = qJDD(4) - t55;
t53 = (-pkin(3) * t72 - pkin(8) * t69) * qJD(2);
t74 = qJD(3) ^ 2;
t75 = qJD(2) ^ 2;
t96 = t65 * t70;
t85 = t73 * t57 + t63 * t96 + t70 * t97;
t26 = -t75 * pkin(2) + qJDD(2) * pkin(7) + t85;
t42 = -t65 * t56 + t67 * t63;
t93 = t72 * t26 + t69 * t42;
t19 = -t74 * pkin(3) + qJDD(3) * pkin(8) + t53 * t90 + t93;
t25 = -qJDD(2) * pkin(2) - t75 * pkin(7) - t101;
t22 = (-t54 - t82) * pkin(8) + (-t55 + t83) * pkin(3) + t25;
t81 = -t68 * t19 + t71 * t22;
t37 = -t61 * mrSges(6,2) + t50 * mrSges(6,3);
t88 = m(6) * (-0.2e1 * qJD(5) * t51 + (t50 * t61 - t32) * qJ(5) + (t50 * t51 + t47) * pkin(4) + t81) + t61 * t37 + t47 * mrSges(6,1);
t10 = m(5) * t81 + t47 * mrSges(5,1) + t61 * t38 + (-t35 - t34) * t51 + (-mrSges(5,3) - mrSges(6,3)) * t32 + t88;
t31 = -t51 * qJD(4) + t71 * qJDD(3) - t68 * t54;
t40 = t61 * mrSges(6,1) - t51 * mrSges(6,3);
t41 = t61 * mrSges(5,1) - t51 * mrSges(5,3);
t39 = t61 * pkin(4) - t51 * qJ(5);
t46 = t50 ^ 2;
t94 = t71 * t19 + t68 * t22;
t87 = m(6) * (-t46 * pkin(4) + t31 * qJ(5) + 0.2e1 * qJD(5) * t50 - t61 * t39 + t94) + t50 * t34 + t31 * mrSges(6,3);
t11 = m(5) * t94 + t31 * mrSges(5,3) + t50 * t35 + (-t41 - t40) * t61 + (-mrSges(5,2) - mrSges(6,2)) * t47 + t87;
t58 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t91;
t59 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t90;
t100 = m(4) * t25 - t55 * mrSges(4,1) + t54 * mrSges(4,2) + t71 * t10 + t68 * t11 + (t69 * t58 - t72 * t59) * qJD(2);
t80 = -t69 * t26 + t72 * t42;
t18 = -qJDD(3) * pkin(3) - t74 * pkin(8) + t53 * t91 - t80;
t86 = m(6) * (-t31 * pkin(4) - t46 * qJ(5) + t51 * t39 + qJDD(5) + t18) + t32 * mrSges(6,2) + t51 * t40;
t99 = m(5) * t18 + t32 * mrSges(5,2) - (mrSges(5,1) + mrSges(6,1)) * t31 + t51 * t41 - (t38 + t37) * t50 + t86;
t8 = m(3) * t101 + qJDD(2) * mrSges(3,1) - t75 * mrSges(3,2) - t100;
t98 = t73 * t8;
t52 = (-mrSges(4,1) * t72 + mrSges(4,2) * t69) * qJD(2);
t12 = m(4) * t80 + qJDD(3) * mrSges(4,1) - t54 * mrSges(4,3) + qJD(3) * t59 - t52 * t91 - t99;
t9 = m(4) * t93 - qJDD(3) * mrSges(4,2) + t55 * mrSges(4,3) - qJD(3) * t58 - t68 * t10 + t71 * t11 + t52 * t90;
t4 = m(3) * t85 - t75 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t69 * t12 + t72 * t9;
t6 = m(3) * t42 + t72 * t12 + t69 * t9;
t84 = m(2) * t63 + t4 * t96 + t67 * t6 + t65 * t98;
t2 = m(2) * t57 + t73 * t4 - t70 * t8;
t1 = m(2) * t56 - t65 * t6 + (t4 * t70 + t98) * t67;
t3 = [-m(1) * g(1) - t64 * t1 + t66 * t2, t2, t4, t9, t11, -t47 * mrSges(6,2) - t61 * t40 + t87; -m(1) * g(2) + t66 * t1 + t64 * t2, t1, t8, t12, t10, -t32 * mrSges(6,3) - t51 * t34 + t88; -m(1) * g(3) + t84, t84, t6, t100, t99, -t31 * mrSges(6,1) - t50 * t37 + t86;];
f_new = t3;
