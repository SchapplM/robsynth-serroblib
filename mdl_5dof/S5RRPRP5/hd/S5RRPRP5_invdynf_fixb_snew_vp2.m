% Calculate vector of cutting forces with Newton-Euler
% S5RRPRP5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRPRP5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP5_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP5_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP5_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP5_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP5_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:54:04
% EndTime: 2019-12-31 19:54:07
% DurationCPUTime: 1.32s
% Computational Cost: add. (13689->161), mult. (31759->210), div. (0->0), fcn. (21142->8), ass. (0->76)
t79 = sin(qJ(2));
t81 = cos(qJ(2));
t99 = qJD(1) * qJD(2);
t65 = qJDD(1) * t79 + t81 * t99;
t66 = qJDD(1) * t81 - t79 * t99;
t101 = qJD(1) * t79;
t68 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t101;
t100 = qJD(1) * t81;
t69 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t100;
t83 = qJD(1) ^ 2;
t76 = sin(pkin(8));
t77 = cos(pkin(8));
t49 = -t65 * t76 + t66 * t77;
t50 = t65 * t77 + t66 * t76;
t59 = (-t76 * t79 + t77 * t81) * qJD(1);
t51 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t59;
t60 = (t76 * t81 + t77 * t79) * qJD(1);
t52 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t60;
t67 = qJD(2) * pkin(2) - qJ(3) * t101;
t75 = t81 ^ 2;
t80 = sin(qJ(1));
t82 = cos(qJ(1));
t96 = t80 * g(1) - t82 * g(2);
t91 = -qJDD(1) * pkin(1) - t96;
t87 = -t66 * pkin(2) + qJDD(3) + t67 * t101 + (-qJ(3) * t75 - pkin(6)) * t83 + t91;
t106 = cos(qJ(4));
t78 = sin(qJ(4));
t43 = t106 * t60 + t59 * t78;
t23 = qJD(4) * t43 - t106 * t49 + t50 * t78;
t42 = -t106 * t59 + t60 * t78;
t24 = -qJD(4) * t42 + t106 * t50 + t49 * t78;
t74 = qJD(2) + qJD(4);
t38 = -mrSges(5,2) * t74 - mrSges(5,3) * t42;
t39 = mrSges(5,1) * t74 - mrSges(5,3) * t43;
t53 = qJD(2) * pkin(3) - pkin(7) * t60;
t58 = t59 ^ 2;
t86 = -t49 * pkin(3) - t58 * pkin(7) + t60 * t53 + t87;
t40 = -mrSges(6,1) * t74 + mrSges(6,2) * t43;
t41 = -mrSges(6,2) * t42 + mrSges(6,3) * t74;
t89 = t24 * mrSges(6,3) + t43 * t40 - m(6) * (-0.2e1 * qJD(5) * t43 + (t42 * t74 - t24) * qJ(5) + (t43 * t74 + t23) * pkin(4) + t86) - t23 * mrSges(6,1) - t42 * t41;
t88 = m(5) * t86 + t23 * mrSges(5,1) + t24 * mrSges(5,2) + t42 * t38 + t43 * t39 - t89;
t85 = -m(4) * t87 + t49 * mrSges(4,1) - t50 * mrSges(4,2) + t59 * t51 - t60 * t52 - t88;
t110 = (t68 * t79 - t69 * t81) * qJD(1) + m(3) * (-t83 * pkin(6) + t91) - t66 * mrSges(3,1) + t65 * mrSges(3,2) - t85;
t93 = -g(1) * t82 - g(2) * t80;
t62 = -pkin(1) * t83 + qJDD(1) * pkin(6) + t93;
t105 = t79 * t62;
t28 = mrSges(6,1) * t42 - mrSges(6,3) * t43;
t102 = -mrSges(5,1) * t42 - mrSges(5,2) * t43 - t28;
t104 = -mrSges(5,3) - mrSges(6,2);
t27 = pkin(4) * t42 - qJ(5) * t43;
t72 = t74 ^ 2;
t73 = qJDD(2) + qJDD(4);
t107 = pkin(2) * t83;
t33 = qJDD(2) * pkin(2) - t65 * qJ(3) - t105 + (qJ(3) * t99 + t107 * t79 - g(3)) * t81;
t95 = -g(3) * t79 + t81 * t62;
t34 = qJ(3) * t66 - qJD(2) * t67 - t107 * t75 + t95;
t94 = -0.2e1 * qJD(3) * t60 + t77 * t33 - t76 * t34;
t16 = (qJD(2) * t59 - t50) * pkin(7) + (t59 * t60 + qJDD(2)) * pkin(3) + t94;
t97 = 0.2e1 * qJD(3) * t59 + t76 * t33 + t77 * t34;
t18 = -pkin(3) * t58 + pkin(7) * t49 - qJD(2) * t53 + t97;
t90 = t106 * t16 - t18 * t78;
t108 = m(6) * (-pkin(4) * t73 - qJ(5) * t72 + t27 * t43 + qJDD(5) - t90);
t10 = m(5) * t90 - t108 + (t38 + t41) * t74 + (mrSges(5,1) + mrSges(6,1)) * t73 + t102 * t43 + t104 * t24;
t46 = -mrSges(4,1) * t59 + mrSges(4,2) * t60;
t103 = t106 * t18 + t16 * t78;
t98 = m(6) * (-pkin(4) * t72 + qJ(5) * t73 + 0.2e1 * qJD(5) * t74 - t27 * t42 + t103) + t74 * t40 + t73 * mrSges(6,3);
t9 = m(5) * t103 - t73 * mrSges(5,2) + t102 * t42 + t104 * t23 - t74 * t39 + t98;
t6 = m(4) * t94 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t50 + qJD(2) * t51 + t10 * t106 - t46 * t60 + t78 * t9;
t64 = (-mrSges(3,1) * t81 + mrSges(3,2) * t79) * qJD(1);
t7 = m(4) * t97 - qJDD(2) * mrSges(4,2) + mrSges(4,3) * t49 - qJD(2) * t52 - t10 * t78 + t106 * t9 + t46 * t59;
t4 = m(3) * (-t81 * g(3) - t105) - t65 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t64 * t101 + qJD(2) * t69 + t76 * t7 + t77 * t6;
t5 = m(3) * t95 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t66 - qJD(2) * t68 + t100 * t64 - t6 * t76 + t7 * t77;
t109 = t4 * t81 + t5 * t79;
t8 = m(2) * t96 + qJDD(1) * mrSges(2,1) - t83 * mrSges(2,2) - t110;
t1 = m(2) * t93 - mrSges(2,1) * t83 - qJDD(1) * mrSges(2,2) - t4 * t79 + t5 * t81;
t2 = [-m(1) * g(1) + t1 * t82 - t8 * t80, t1, t5, t7, t9, -t23 * mrSges(6,2) - t42 * t28 + t98; -m(1) * g(2) + t1 * t80 + t8 * t82, t8, t4, t6, t10, -t89; (-m(1) - m(2)) * g(3) + t109, -m(2) * g(3) + t109, t110, -t85, t88, -t73 * mrSges(6,1) + t24 * mrSges(6,2) + t43 * t28 - t74 * t41 + t108;];
f_new = t2;
