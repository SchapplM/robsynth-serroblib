% Calculate vector of cutting forces with Newton-Euler
% S5RPRPR1
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-05 17:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRPR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR1_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR1_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR1_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR1_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR1_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:47:06
% EndTime: 2019-12-05 17:47:08
% DurationCPUTime: 0.97s
% Computational Cost: add. (9516->137), mult. (20663->181), div. (0->0), fcn. (12967->8), ass. (0->71)
t71 = sin(qJ(1));
t74 = cos(qJ(1));
t85 = -t74 * g(1) - t71 * g(2);
t83 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t85;
t100 = -m(2) - m(3);
t99 = -pkin(1) - pkin(6);
t98 = (mrSges(2,1) - mrSges(3,2));
t97 = -mrSges(2,2) + mrSges(3,3);
t75 = qJD(1) ^ 2;
t89 = t71 * g(1) - t74 * g(2);
t81 = -t75 * qJ(2) + qJDD(2) - t89;
t42 = t99 * qJDD(1) + t81;
t70 = sin(qJ(3));
t73 = cos(qJ(3));
t96 = t70 * g(3) + t73 * t42;
t95 = qJD(1) * t70;
t94 = qJD(1) * t73;
t93 = qJD(1) * qJD(3);
t87 = t70 * t93;
t56 = t73 * qJDD(1) - t87;
t22 = (-t56 - t87) * qJ(4) + (-t70 * t73 * t75 + qJDD(3)) * pkin(3) + t96;
t55 = -t70 * qJDD(1) - t73 * t93;
t58 = qJD(3) * pkin(3) - qJ(4) * t94;
t66 = t70 ^ 2;
t88 = -t73 * g(3) + t70 * t42;
t23 = -t66 * t75 * pkin(3) + t55 * qJ(4) - qJD(3) * t58 + t88;
t67 = sin(pkin(8));
t68 = cos(pkin(8));
t48 = (-t67 * t73 - t68 * t70) * qJD(1);
t91 = 0.2e1 * qJD(4) * t48 + t67 * t22 + t68 * t23;
t49 = (-t67 * t70 + t68 * t73) * qJD(1);
t31 = -t48 * mrSges(5,1) + t49 * mrSges(5,2);
t34 = t67 * t55 + t68 * t56;
t39 = -qJD(3) * mrSges(5,2) + t48 * mrSges(5,3);
t69 = sin(qJ(5));
t72 = cos(qJ(5));
t86 = -0.2e1 * qJD(4) * t49 + t68 * t22 - t67 * t23;
t10 = (qJD(3) * t48 - t34) * pkin(7) + (t48 * t49 + qJDD(3)) * pkin(4) + t86;
t33 = t68 * t55 - t67 * t56;
t41 = qJD(3) * pkin(4) - t49 * pkin(7);
t47 = t48 ^ 2;
t11 = -t47 * pkin(4) + t33 * pkin(7) - qJD(3) * t41 + t91;
t28 = t72 * t48 - t69 * t49;
t16 = t28 * qJD(5) + t69 * t33 + t72 * t34;
t29 = t69 * t48 + t72 * t49;
t18 = -t28 * mrSges(6,1) + t29 * mrSges(6,2);
t64 = qJD(3) + qJD(5);
t26 = -t64 * mrSges(6,2) + t28 * mrSges(6,3);
t63 = qJDD(3) + qJDD(5);
t8 = m(6) * (t72 * t10 - t69 * t11) - t16 * mrSges(6,3) + t63 * mrSges(6,1) - t29 * t18 + t64 * t26;
t15 = -t29 * qJD(5) + t72 * t33 - t69 * t34;
t27 = t64 * mrSges(6,1) - t29 * mrSges(6,3);
t9 = m(6) * (t69 * t10 + t72 * t11) + t15 * mrSges(6,3) - t63 * mrSges(6,2) + t28 * t18 - t64 * t27;
t5 = m(5) * t86 + qJDD(3) * mrSges(5,1) - t34 * mrSges(5,3) + qJD(3) * t39 - t49 * t31 + t69 * t9 + t72 * t8;
t54 = (mrSges(4,1) * t70 + mrSges(4,2) * t73) * qJD(1);
t57 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t95;
t40 = qJD(3) * mrSges(5,1) - t49 * mrSges(5,3);
t6 = m(5) * t91 - qJDD(3) * mrSges(5,2) + t33 * mrSges(5,3) - qJD(3) * t40 + t48 * t31 - t69 * t8 + t72 * t9;
t3 = m(4) * t96 + qJDD(3) * mrSges(4,1) - t56 * mrSges(4,3) + qJD(3) * t57 + t68 * t5 - t54 * t94 + t67 * t6;
t59 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t94;
t4 = m(4) * t88 - qJDD(3) * mrSges(4,2) + t55 * mrSges(4,3) - qJD(3) * t59 - t67 * t5 - t54 * t95 + t68 * t6;
t90 = -t70 * t3 + t73 * t4;
t82 = -m(3) * (-qJDD(1) * pkin(1) + t81) - t73 * t3 - t70 * t4;
t78 = -t55 * pkin(3) + qJDD(4) + t58 * t94 + (-qJ(4) * t66 + t99) * t75 + t83;
t80 = -t15 * mrSges(6,1) - t28 * t26 + m(6) * (-t33 * pkin(4) - t47 * pkin(7) + t49 * t41 + t78) + t16 * mrSges(6,2) + t29 * t27;
t79 = m(5) * t78 - t33 * mrSges(5,1) + t34 * mrSges(5,2) - t48 * t39 + t49 * t40 + t80;
t77 = -t55 * mrSges(4,1) + m(4) * (t99 * t75 + t83) + t57 * t95 + t59 * t94 + t56 * mrSges(4,2) + t79;
t76 = -m(3) * (t75 * pkin(1) - t83) + t77;
t7 = m(2) * t85 + t97 * qJDD(1) - (t98 * t75) + t76;
t1 = m(2) * t89 + t98 * qJDD(1) + t97 * t75 + t82;
t2 = [-m(1) * g(1) - t71 * t1 + t74 * t7, t7, -m(3) * g(3) + t90, t4, t6, t9; -m(1) * g(2) + t74 * t1 + t71 * t7, t1, -(t75 * mrSges(3,2)) - qJDD(1) * mrSges(3,3) - t76, t3, t5, t8; (-m(1) + t100) * g(3) + t90, t100 * g(3) + t90, qJDD(1) * mrSges(3,2) - t75 * mrSges(3,3) - t82, t77, t79, t80;];
f_new = t2;
