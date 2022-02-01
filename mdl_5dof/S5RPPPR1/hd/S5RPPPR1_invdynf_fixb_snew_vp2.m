% Calculate vector of cutting forces with Newton-Euler
% S5RPPPR1
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% m [6x1]
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
% Datum: 2022-01-20 09:13
% Revision: 008671b0a00594318b890887636eaaff83cd5e2f (2021-12-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPPPR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR1_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR1_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR1_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR1_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:12:19
% EndTime: 2022-01-20 09:12:21
% DurationCPUTime: 0.98s
% Computational Cost: add. (8557->119), mult. (19883->176), div. (0->0), fcn. (12023->10), ass. (0->78)
t71 = qJD(1) ^ 2;
t68 = sin(qJ(1));
t70 = cos(qJ(1));
t89 = t68 * g(1) - t70 * g(2);
t46 = qJDD(1) * pkin(1) + t89;
t85 = -t70 * g(1) - t68 * g(2);
t47 = -t71 * pkin(1) + t85;
t63 = sin(pkin(7));
t66 = cos(pkin(7));
t87 = t66 * t46 - t63 * t47;
t74 = -t71 * qJ(3) + qJDD(3) - t87;
t62 = sin(pkin(8));
t65 = cos(pkin(8));
t82 = -pkin(3) * t65 - qJ(4) * t62;
t97 = qJD(1) * t62;
t105 = (-pkin(2) + t82) * qJDD(1) + t74 - 0.2e1 * qJD(4) * t97;
t58 = t62 ^ 2;
t104 = (t65 ^ 2 + t58) * mrSges(4,3);
t103 = pkin(6) * t62;
t61 = sin(pkin(9));
t102 = t61 * t62;
t64 = cos(pkin(9));
t101 = t62 * t64;
t100 = t105 * t64;
t99 = t63 * t46 + t66 * t47;
t96 = t65 * qJD(1);
t95 = qJDD(1) * t61;
t84 = -t65 * mrSges(4,1) + t62 * mrSges(4,2);
t45 = t84 * qJD(1);
t44 = t82 * qJD(1);
t31 = -t71 * pkin(2) + qJDD(1) * qJ(3) + t99;
t60 = -g(3) + qJDD(2);
t88 = -t62 * t31 + t65 * t60;
t90 = 0.2e1 * qJD(1) * qJD(3);
t18 = t44 * t97 + t62 * t90 + qJDD(4) - t88;
t67 = sin(qJ(5));
t69 = cos(qJ(5));
t76 = (-t61 * t67 + t64 * t69) * t62;
t35 = qJD(1) * t76;
t77 = (-t61 * t69 - t64 * t67) * t62;
t26 = -t35 * qJD(5) + qJDD(1) * t77;
t34 = qJD(1) * t77;
t27 = t34 * qJD(5) + qJDD(1) * t76;
t51 = qJD(5) - t96;
t32 = -t51 * mrSges(6,2) + t34 * mrSges(6,3);
t33 = t51 * mrSges(6,1) - t35 * mrSges(6,3);
t80 = -pkin(4) * t65 - pkin(6) * t101;
t43 = t80 * qJD(1);
t93 = t61 ^ 2 * t58 * t71;
t73 = t26 * mrSges(6,1) + t34 * t32 - m(6) * (-pkin(6) * t93 + (qJD(1) * t43 * t64 + pkin(4) * t95) * t62 + t18) - t35 * t33 - t27 * mrSges(6,2);
t72 = m(5) * t18 - t73;
t78 = t65 * mrSges(5,2) - mrSges(5,3) * t102;
t38 = t78 * qJD(1);
t79 = -t65 * mrSges(5,1) - mrSges(5,3) * t101;
t39 = t79 * qJD(1);
t81 = t38 * t61 + t39 * t64;
t83 = mrSges(5,1) * t61 + mrSges(5,2) * t64;
t10 = m(4) * t88 + ((-mrSges(4,3) - t83) * qJDD(1) + (-0.2e1 * m(4) * qJD(3) - t45 - t81) * qJD(1)) * t62 - t72;
t91 = t62 * t60 + (t31 + t90) * t65;
t19 = t44 * t96 + t91;
t13 = t80 * qJDD(1) + (-t19 + (-pkin(4) * t58 * t64 + t65 * t103) * t71) * t61 + t100;
t92 = t105 * t61 + t64 * t19;
t14 = -pkin(4) * t93 - t95 * t103 + t43 * t96 + t92;
t24 = -t34 * mrSges(6,1) + t35 * mrSges(6,2);
t50 = -t65 * qJDD(1) + qJDD(5);
t11 = m(6) * (t69 * t13 - t67 * t14) - t27 * mrSges(6,3) + t50 * mrSges(6,1) - t35 * t24 + t51 * t32;
t12 = m(6) * (t67 * t13 + t69 * t14) + t26 * mrSges(6,3) - t50 * mrSges(6,2) + t34 * t24 - t51 * t33;
t36 = t83 * t97;
t7 = m(5) * (-t61 * t19 + t100) + t67 * t12 + t69 * t11 + t79 * qJDD(1) + (-t36 * t101 - t65 * t38) * qJD(1);
t8 = m(5) * t92 + t69 * t12 - t67 * t11 + t78 * qJDD(1) + (-t36 * t102 + t65 * t39) * qJD(1);
t6 = m(4) * t91 + t64 * t8 - t61 * t7 + (qJDD(1) * mrSges(4,3) + qJD(1) * t45) * t65;
t94 = m(3) * t60 + t65 * t10 + t62 * t6;
t75 = m(4) * (-qJDD(1) * pkin(2) + t74) + t61 * t8 + t64 * t7;
t4 = m(3) * t87 + (-mrSges(3,2) + t104) * t71 + (mrSges(3,1) - t84) * qJDD(1) - t75;
t3 = m(3) * t99 - t71 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t62 * t10 + t65 * t6;
t2 = m(2) * t85 - t71 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t66 * t3 - t63 * t4;
t1 = m(2) * t89 + qJDD(1) * mrSges(2,1) - t71 * mrSges(2,2) + t63 * t3 + t66 * t4;
t5 = [-m(1) * g(1) - t68 * t1 + t70 * t2, t2, t3, t6, t8, t12; -m(1) * g(2) + t70 * t1 + t68 * t2, t1, t4, t10, t7, t11; (-m(1) - m(2)) * g(3) + t94, -m(2) * g(3) + t94, t94, t84 * qJDD(1) - t71 * t104 + t75, (t81 * qJD(1) + t83 * qJDD(1)) * t62 + t72, -t73;];
f_new = t5;
