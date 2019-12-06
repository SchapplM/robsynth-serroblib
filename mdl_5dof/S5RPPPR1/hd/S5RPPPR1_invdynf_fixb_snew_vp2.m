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
% Datum: 2019-12-05 17:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:28:50
% EndTime: 2019-12-05 17:28:51
% DurationCPUTime: 0.90s
% Computational Cost: add. (8557->119), mult. (19883->176), div. (0->0), fcn. (12023->10), ass. (0->78)
t73 = qJD(1) ^ 2;
t70 = sin(qJ(1));
t72 = cos(qJ(1));
t100 = t72 * g(2) + t70 * g(3);
t46 = qJDD(1) * pkin(1) + t100;
t90 = t70 * g(2) - t72 * g(3);
t47 = -t73 * pkin(1) + t90;
t65 = sin(pkin(7));
t68 = cos(pkin(7));
t88 = t68 * t46 - t65 * t47;
t76 = -t73 * qJ(3) + qJDD(3) - t88;
t64 = sin(pkin(8));
t67 = cos(pkin(8));
t84 = -pkin(3) * t67 - qJ(4) * t64;
t98 = qJD(1) * t64;
t107 = (-pkin(2) + t84) * qJDD(1) + t76 - 0.2e1 * qJD(4) * t98;
t60 = t64 ^ 2;
t106 = (t67 ^ 2 + t60) * mrSges(4,3);
t105 = pkin(6) * t64;
t63 = sin(pkin(9));
t104 = t63 * t64;
t66 = cos(pkin(9));
t103 = t64 * t66;
t102 = t107 * t66;
t101 = t65 * t46 + t68 * t47;
t97 = t67 * qJD(1);
t96 = qJDD(1) * t63;
t86 = -t67 * mrSges(4,1) + t64 * mrSges(4,2);
t45 = t86 * qJD(1);
t44 = t84 * qJD(1);
t31 = -t73 * pkin(2) + qJDD(1) * qJ(3) + t101;
t62 = -g(1) + qJDD(2);
t89 = -t64 * t31 + t67 * t62;
t91 = 0.2e1 * qJD(1) * qJD(3);
t18 = t44 * t98 + t64 * t91 + qJDD(4) - t89;
t69 = sin(qJ(5));
t71 = cos(qJ(5));
t78 = (-t63 * t69 + t66 * t71) * t64;
t35 = qJD(1) * t78;
t79 = (-t63 * t71 - t66 * t69) * t64;
t26 = -t35 * qJD(5) + qJDD(1) * t79;
t34 = qJD(1) * t79;
t27 = t34 * qJD(5) + qJDD(1) * t78;
t51 = qJD(5) - t97;
t32 = -t51 * mrSges(6,2) + t34 * mrSges(6,3);
t33 = t51 * mrSges(6,1) - t35 * mrSges(6,3);
t82 = -pkin(4) * t67 - pkin(6) * t103;
t43 = t82 * qJD(1);
t94 = t63 ^ 2 * t60 * t73;
t75 = t26 * mrSges(6,1) + t34 * t32 - m(6) * (-pkin(6) * t94 + (qJD(1) * t43 * t66 + pkin(4) * t96) * t64 + t18) - t35 * t33 - t27 * mrSges(6,2);
t74 = m(5) * t18 - t75;
t80 = t67 * mrSges(5,2) - mrSges(5,3) * t104;
t38 = t80 * qJD(1);
t81 = -t67 * mrSges(5,1) - mrSges(5,3) * t103;
t39 = t81 * qJD(1);
t83 = t38 * t63 + t39 * t66;
t85 = mrSges(5,1) * t63 + mrSges(5,2) * t66;
t10 = m(4) * t89 + ((-mrSges(4,3) - t85) * qJDD(1) + (-0.2e1 * m(4) * qJD(3) - t45 - t83) * qJD(1)) * t64 - t74;
t92 = t64 * t62 + (t31 + t91) * t67;
t19 = t44 * t97 + t92;
t13 = t82 * qJDD(1) + (-t19 + (-pkin(4) * t60 * t66 + t67 * t105) * t73) * t63 + t102;
t93 = t107 * t63 + t66 * t19;
t14 = -pkin(4) * t94 - t96 * t105 + t43 * t97 + t93;
t24 = -t34 * mrSges(6,1) + t35 * mrSges(6,2);
t50 = -t67 * qJDD(1) + qJDD(5);
t11 = m(6) * (t71 * t13 - t69 * t14) - t27 * mrSges(6,3) + t50 * mrSges(6,1) - t35 * t24 + t51 * t32;
t12 = m(6) * (t69 * t13 + t71 * t14) + t26 * mrSges(6,3) - t50 * mrSges(6,2) + t34 * t24 - t51 * t33;
t36 = t85 * t98;
t7 = m(5) * (-t63 * t19 + t102) + t69 * t12 + t71 * t11 + t81 * qJDD(1) + (-t36 * t103 - t67 * t38) * qJD(1);
t8 = m(5) * t93 + t71 * t12 - t69 * t11 + t80 * qJDD(1) + (-t36 * t104 + t67 * t39) * qJD(1);
t6 = m(4) * t92 + t66 * t8 - t63 * t7 + (qJDD(1) * mrSges(4,3) + qJD(1) * t45) * t67;
t95 = m(3) * t62 + t67 * t10 + t64 * t6;
t77 = m(4) * (-qJDD(1) * pkin(2) + t76) + t63 * t8 + t66 * t7;
t4 = m(3) * t88 + (-mrSges(3,2) + t106) * t73 + (mrSges(3,1) - t86) * qJDD(1) - t77;
t3 = m(3) * t101 - t73 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t64 * t10 + t67 * t6;
t2 = m(2) * t100 + qJDD(1) * mrSges(2,1) - t73 * mrSges(2,2) + t65 * t3 + t68 * t4;
t1 = m(2) * t90 - t73 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t68 * t3 - t65 * t4;
t5 = [(-m(1) - m(2)) * g(1) + t95, t1, t3, t6, t8, t12; -m(1) * g(2) - t70 * t1 - t72 * t2, t2, t4, t10, t7, t11; -m(1) * g(3) + t72 * t1 - t70 * t2, -m(2) * g(1) + t95, t95, t86 * qJDD(1) - t73 * t106 + t77, (t83 * qJD(1) + t85 * qJDD(1)) * t64 + t74, -t75;];
f_new = t5;
