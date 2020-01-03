% Calculate vector of cutting forces with Newton-Euler
% S5RPPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2020-01-03 11:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPPRP1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP1_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP1_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP1_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:25:12
% EndTime: 2020-01-03 11:25:15
% DurationCPUTime: 0.60s
% Computational Cost: add. (4331->121), mult. (9010->161), div. (0->0), fcn. (4988->8), ass. (0->70)
t59 = sin(pkin(8));
t55 = t59 ^ 2;
t61 = cos(pkin(8));
t101 = (t61 ^ 2 + t55) * mrSges(4,3);
t74 = -pkin(3) * t61 - pkin(6) * t59;
t44 = t74 * qJD(1);
t67 = qJD(1) ^ 2;
t64 = sin(qJ(1));
t66 = cos(qJ(1));
t73 = -t66 * g(2) - t64 * g(3);
t45 = qJDD(1) * pkin(1) + t73;
t80 = -t64 * g(2) + t66 * g(3);
t46 = -t67 * pkin(1) + t80;
t60 = sin(pkin(7));
t62 = cos(pkin(7));
t95 = t60 * t45 + t62 * t46;
t24 = -t67 * pkin(2) + qJDD(1) * qJ(3) + t95;
t58 = -g(1) + qJDD(2);
t79 = -t59 * t24 + t61 * t58;
t82 = 0.2e1 * qJD(1) * qJD(3);
t93 = qJD(1) * t59;
t16 = t44 * t93 + t59 * t82 - t79;
t63 = sin(qJ(4));
t65 = cos(qJ(4));
t90 = qJD(1) * qJD(4);
t37 = (-qJDD(1) * t63 - t65 * t90) * t59;
t38 = (qJDD(1) * t65 - t63 * t90) * t59;
t92 = t61 * qJD(1);
t49 = qJD(4) - t92;
t85 = t63 * t93;
t30 = -t49 * mrSges(6,2) - mrSges(6,3) * t85;
t81 = qJ(5) * t93;
t32 = t49 * pkin(4) - t65 * t81;
t84 = t65 * t93;
t33 = t49 * mrSges(6,1) - mrSges(6,3) * t84;
t99 = t55 * t67;
t88 = t63 ^ 2 * t99;
t77 = m(6) * (-t37 * pkin(4) - qJ(5) * t88 + t32 * t84 + qJDD(5) + t16) + t30 * t85 + t33 * t84 + t38 * mrSges(6,2);
t100 = -m(5) * t16 - t38 * mrSges(5,2) + (mrSges(5,1) + mrSges(6,1)) * t37 - t77;
t76 = -0.2e1 * qJD(5) * t93;
t86 = t59 * t58 + (t24 + t82) * t61;
t17 = t44 * t92 + t86;
t78 = t62 * t45 - t60 * t46;
t69 = -t67 * qJ(3) + qJDD(3) - t78;
t20 = (-pkin(2) + t74) * qJDD(1) + t69;
t96 = t65 * t17 + t63 * t20;
t97 = m(6) * (-pkin(4) * t88 + t37 * qJ(5) - t49 * t32 + t63 * t76 + t96) + t37 * mrSges(6,3);
t91 = qJDD(1) * mrSges(4,3);
t72 = -t61 * mrSges(4,1) + t59 * mrSges(4,2);
t43 = t72 * qJD(1);
t31 = -t49 * mrSges(5,2) - mrSges(5,3) * t85;
t34 = t49 * mrSges(5,1) - mrSges(5,3) * t84;
t71 = t31 * t63 + t34 * t65;
t10 = m(4) * t79 + (-t91 + (-0.2e1 * m(4) * qJD(3) - t43 - t71) * qJD(1)) * t59 + t100;
t19 = t65 * t20;
t48 = -t61 * qJDD(1) + qJDD(4);
t35 = (mrSges(6,1) * t63 + mrSges(6,2) * t65) * t93;
t75 = (-t35 - (mrSges(5,1) * t63 + mrSges(5,2) * t65) * t93) * t93;
t87 = m(6) * (t65 * t76 + t48 * pkin(4) - t38 * qJ(5) + t19 + (-pkin(4) * t65 * t99 - t49 * t81 - t17) * t63) + t49 * t30 + t48 * mrSges(6,1);
t7 = m(5) * (-t63 * t17 + t19) + t48 * mrSges(5,1) + t49 * t31 + (-mrSges(5,3) - mrSges(6,3)) * t38 + t65 * t75 + t87;
t8 = m(5) * t96 + t37 * mrSges(5,3) + (-t34 - t33) * t49 + (-mrSges(5,2) - mrSges(6,2)) * t48 + t63 * t75 + t97;
t6 = m(4) * t86 + t65 * t8 - t63 * t7 + (qJD(1) * t43 + t91) * t61;
t89 = m(3) * t58 + t61 * t10 + t59 * t6;
t83 = t35 * t93;
t70 = m(4) * (-qJDD(1) * pkin(2) + t69) + t63 * t8 + t65 * t7;
t4 = m(3) * t78 + (-mrSges(3,2) + t101) * t67 + (mrSges(3,1) - t72) * qJDD(1) - t70;
t3 = m(3) * t95 - t67 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t59 * t10 + t61 * t6;
t2 = m(2) * t73 + qJDD(1) * mrSges(2,1) - t67 * mrSges(2,2) + t60 * t3 + t62 * t4;
t1 = m(2) * t80 - t67 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t62 * t3 - t60 * t4;
t5 = [(-m(1) - m(2)) * g(1) + t89, t1, t3, t6, t8, -t48 * mrSges(6,2) - t49 * t33 - t63 * t83 + t97; -m(1) * g(2) + t64 * t1 + t66 * t2, t2, t4, t10, t7, -t38 * mrSges(6,3) - t65 * t83 + t87; -m(1) * g(3) - t66 * t1 + t64 * t2, -m(2) * g(1) + t89, t89, t72 * qJDD(1) - t67 * t101 + t70, t71 * t93 - t100, -t37 * mrSges(6,1) + t77;];
f_new = t5;
