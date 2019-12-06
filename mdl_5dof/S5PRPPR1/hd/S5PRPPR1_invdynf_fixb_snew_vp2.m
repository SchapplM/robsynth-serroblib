% Calculate vector of cutting forces with Newton-Euler
% S5PRPPR1
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
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
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRPPR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR1_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR1_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR1_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR1_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR1_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:21:52
% EndTime: 2019-12-05 15:21:54
% DurationCPUTime: 0.91s
% Computational Cost: add. (7593->112), mult. (18434->170), div. (0->0), fcn. (12023->10), ass. (0->77)
t71 = qJD(2) ^ 2;
t63 = sin(pkin(7));
t66 = cos(pkin(7));
t48 = t63 * g(1) - t66 * g(2);
t49 = -t66 * g(1) - t63 * g(2);
t68 = sin(qJ(2));
t70 = cos(qJ(2));
t86 = t70 * t48 - t68 * t49;
t74 = -t71 * qJ(3) + qJDD(3) - t86;
t62 = sin(pkin(8));
t65 = cos(pkin(8));
t82 = -pkin(3) * t65 - qJ(4) * t62;
t96 = qJD(2) * t62;
t104 = (-pkin(2) + t82) * qJDD(2) + t74 - 0.2e1 * qJD(4) * t96;
t58 = t62 ^ 2;
t103 = (t65 ^ 2 + t58) * mrSges(4,3);
t102 = pkin(6) * t62;
t61 = sin(pkin(9));
t101 = t61 * t62;
t64 = cos(pkin(9));
t100 = t62 * t64;
t99 = t104 * t64;
t98 = t68 * t48 + t70 * t49;
t95 = t65 * qJD(2);
t94 = qJDD(2) * t61;
t84 = -t65 * mrSges(4,1) + t62 * mrSges(4,2);
t45 = t84 * qJD(2);
t44 = t82 * qJD(2);
t31 = -t71 * pkin(2) + qJDD(2) * qJ(3) + t98;
t60 = -g(3) + qJDD(1);
t87 = -t62 * t31 + t65 * t60;
t89 = 0.2e1 * qJD(2) * qJD(3);
t18 = t44 * t96 + t62 * t89 + qJDD(4) - t87;
t67 = sin(qJ(5));
t69 = cos(qJ(5));
t76 = (-t61 * t67 + t64 * t69) * t62;
t35 = qJD(2) * t76;
t77 = (-t61 * t69 - t64 * t67) * t62;
t26 = -t35 * qJD(5) + qJDD(2) * t77;
t34 = qJD(2) * t77;
t27 = t34 * qJD(5) + qJDD(2) * t76;
t51 = qJD(5) - t95;
t32 = -t51 * mrSges(6,2) + t34 * mrSges(6,3);
t33 = t51 * mrSges(6,1) - t35 * mrSges(6,3);
t80 = -pkin(4) * t65 - pkin(6) * t100;
t40 = t80 * qJD(2);
t92 = t61 ^ 2 * t58 * t71;
t73 = t26 * mrSges(6,1) + t34 * t32 - m(6) * (-pkin(6) * t92 + (qJD(2) * t40 * t64 + pkin(4) * t94) * t62 + t18) - t35 * t33 - t27 * mrSges(6,2);
t72 = m(5) * t18 - t73;
t78 = t65 * mrSges(5,2) - mrSges(5,3) * t101;
t38 = t78 * qJD(2);
t79 = -t65 * mrSges(5,1) - mrSges(5,3) * t100;
t39 = t79 * qJD(2);
t81 = t38 * t61 + t39 * t64;
t83 = mrSges(5,1) * t61 + mrSges(5,2) * t64;
t10 = m(4) * t87 + ((-mrSges(4,3) - t83) * qJDD(2) + (-0.2e1 * m(4) * qJD(3) - t45 - t81) * qJD(2)) * t62 - t72;
t90 = t62 * t60 + (t31 + t89) * t65;
t19 = t44 * t95 + t90;
t13 = t80 * qJDD(2) + (-t19 + (-pkin(4) * t58 * t64 + t65 * t102) * t71) * t61 + t99;
t91 = t104 * t61 + t64 * t19;
t14 = -pkin(4) * t92 - t94 * t102 + t40 * t95 + t91;
t23 = -t34 * mrSges(6,1) + t35 * mrSges(6,2);
t50 = -t65 * qJDD(2) + qJDD(5);
t11 = m(6) * (t69 * t13 - t67 * t14) - t27 * mrSges(6,3) + t50 * mrSges(6,1) - t35 * t23 + t51 * t32;
t12 = m(6) * (t67 * t13 + t69 * t14) + t26 * mrSges(6,3) - t50 * mrSges(6,2) + t34 * t23 - t51 * t33;
t36 = t83 * t96;
t7 = m(5) * (-t61 * t19 + t99) + t67 * t12 + t69 * t11 + t79 * qJDD(2) + (-t36 * t100 - t65 * t38) * qJD(2);
t8 = m(5) * t91 + t69 * t12 - t67 * t11 + t78 * qJDD(2) + (-t36 * t101 + t65 * t39) * qJD(2);
t6 = m(4) * t90 + t64 * t8 - t61 * t7 + (qJDD(2) * mrSges(4,3) + qJD(2) * t45) * t65;
t93 = m(3) * t60 + t65 * t10 + t62 * t6;
t88 = m(2) * t60 + t93;
t75 = m(4) * (-qJDD(2) * pkin(2) + t74) + t61 * t8 + t64 * t7;
t4 = m(3) * t86 + (-mrSges(3,2) + t103) * t71 + (mrSges(3,1) - t84) * qJDD(2) - t75;
t3 = m(3) * t98 - t71 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t62 * t10 + t65 * t6;
t2 = m(2) * t49 + t70 * t3 - t68 * t4;
t1 = m(2) * t48 + t68 * t3 + t70 * t4;
t5 = [-m(1) * g(1) - t63 * t1 + t66 * t2, t2, t3, t6, t8, t12; -m(1) * g(2) + t66 * t1 + t63 * t2, t1, t4, t10, t7, t11; -m(1) * g(3) + t88, t88, t93, t84 * qJDD(2) - t103 * t71 + t75, (t81 * qJD(2) + t83 * qJDD(2)) * t62 + t72, -t73;];
f_new = t5;
