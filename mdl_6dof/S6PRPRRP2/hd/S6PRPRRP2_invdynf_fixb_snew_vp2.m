% Calculate vector of cutting forces with Newton-Euler
% S6PRPRRP2
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% f_new [3x7]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 23:40
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRPRRP2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP2_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP2_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP2_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP2_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP2_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:35:01
% EndTime: 2019-05-04 23:35:04
% DurationCPUTime: 1.52s
% Computational Cost: add. (19625->146), mult. (35705->190), div. (0->0), fcn. (24210->12), ass. (0->82)
t116 = cos(qJ(5));
t84 = cos(qJ(4));
t104 = t84 * qJD(2);
t77 = sin(pkin(6));
t85 = cos(qJ(2));
t113 = t77 * t85;
t76 = sin(pkin(10));
t79 = cos(pkin(10));
t62 = t76 * g(1) - t79 * g(2);
t80 = cos(pkin(6));
t115 = t62 * t80;
t63 = -t79 * g(1) - t76 * g(2);
t74 = -g(3) + qJDD(1);
t83 = sin(qJ(2));
t92 = t74 * t113 + t85 * t115 - t83 * t63;
t33 = qJDD(2) * pkin(2) + t92;
t87 = qJD(2) ^ 2;
t114 = t77 * t83;
t99 = t74 * t114 + t83 * t115 + t85 * t63;
t34 = -t87 * pkin(2) + t99;
t75 = sin(pkin(11));
t78 = cos(pkin(11));
t108 = t75 * t33 + t78 * t34;
t29 = -t87 * pkin(3) + qJDD(2) * pkin(8) + t108;
t93 = -t77 * t62 + t80 * t74;
t49 = qJDD(3) + t93;
t82 = sin(qJ(4));
t109 = t84 * t29 + t82 * t49;
t59 = (-pkin(4) * t84 - pkin(9) * t82) * qJD(2);
t86 = qJD(4) ^ 2;
t23 = -t86 * pkin(4) + qJDD(4) * pkin(9) + t104 * t59 + t109;
t94 = t78 * t33 - t75 * t34;
t28 = -qJDD(2) * pkin(3) - t87 * pkin(8) - t94;
t103 = qJD(2) * qJD(4);
t96 = t84 * t103;
t60 = t82 * qJDD(2) + t96;
t97 = t82 * t103;
t61 = t84 * qJDD(2) - t97;
t25 = (-t60 - t96) * pkin(9) + (-t61 + t97) * pkin(4) + t28;
t81 = sin(qJ(5));
t110 = t116 * t23 + t81 * t25;
t105 = qJD(2) * t82;
t56 = -qJD(4) * t116 + t105 * t81;
t57 = t81 * qJD(4) + t105 * t116;
t40 = t56 * pkin(5) - t57 * qJ(6);
t69 = qJD(5) - t104;
t47 = -t69 * mrSges(7,1) + t57 * mrSges(7,2);
t54 = qJDD(5) - t61;
t68 = t69 ^ 2;
t102 = m(7) * (-t68 * pkin(5) + t54 * qJ(6) + 0.2e1 * qJD(6) * t69 - t56 * t40 + t110) + t69 * t47 + t54 * mrSges(7,3);
t41 = t56 * mrSges(7,1) - t57 * mrSges(7,3);
t107 = -t56 * mrSges(6,1) - t57 * mrSges(6,2) - t41;
t111 = -mrSges(6,3) - mrSges(7,2);
t37 = t57 * qJD(5) - qJDD(4) * t116 + t81 * t60;
t46 = t69 * mrSges(6,1) - t57 * mrSges(6,3);
t13 = m(6) * t110 - t54 * mrSges(6,2) + t107 * t56 + t111 * t37 - t69 * t46 + t102;
t90 = t116 * t25 - t81 * t23;
t117 = m(7) * (-t54 * pkin(5) - t68 * qJ(6) + t57 * t40 + qJDD(6) - t90);
t38 = -t56 * qJD(5) + t81 * qJDD(4) + t116 * t60;
t45 = -t69 * mrSges(6,2) - t56 * mrSges(6,3);
t48 = -t56 * mrSges(7,2) + t69 * mrSges(7,3);
t14 = m(6) * t90 - t117 + (t45 + t48) * t69 + t107 * t57 + (mrSges(6,1) + mrSges(7,1)) * t54 + t111 * t38;
t64 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t105;
t65 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t104;
t119 = m(5) * t28 - t61 * mrSges(5,1) + t60 * mrSges(5,2) + (t64 * t82 - t65 * t84) * qJD(2) + t116 * t14 + t81 * t13;
t95 = -t82 * t29 + t84 * t49;
t22 = -qJDD(4) * pkin(4) - t86 * pkin(9) + t59 * t105 - t95;
t101 = m(7) * (-0.2e1 * qJD(6) * t57 + (t56 * t69 - t38) * qJ(6) + (t57 * t69 + t37) * pkin(5) + t22) + t37 * mrSges(7,1) + t56 * t48;
t118 = m(6) * t22 + t37 * mrSges(6,1) + (t46 - t47) * t57 + (mrSges(6,2) - mrSges(7,3)) * t38 + t56 * t45 + t101;
t58 = (-mrSges(5,1) * t84 + mrSges(5,2) * t82) * qJD(2);
t12 = m(5) * t109 - qJDD(4) * mrSges(5,2) + t61 * mrSges(5,3) - qJD(4) * t64 + t104 * t58 + t116 * t13 - t81 * t14;
t16 = m(5) * t95 + qJDD(4) * mrSges(5,1) - t60 * mrSges(5,3) + qJD(4) * t65 - t105 * t58 - t118;
t100 = m(4) * t49 + t82 * t12 + t84 * t16;
t10 = m(4) * t94 + qJDD(2) * mrSges(4,1) - t87 * mrSges(4,2) - t119;
t7 = m(4) * t108 - t87 * mrSges(4,1) - qJDD(2) * mrSges(4,2) + t84 * t12 - t82 * t16;
t5 = m(3) * t92 + qJDD(2) * mrSges(3,1) - t87 * mrSges(3,2) + t78 * t10 + t75 * t7;
t6 = m(3) * t99 - t87 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t75 * t10 + t78 * t7;
t9 = m(3) * t93 + t100;
t98 = m(2) * t74 + t5 * t113 + t6 * t114 + t80 * t9;
t2 = m(2) * t63 - t83 * t5 + t85 * t6;
t1 = m(2) * t62 - t77 * t9 + (t5 * t85 + t6 * t83) * t80;
t3 = [-m(1) * g(1) - t76 * t1 + t79 * t2, t2, t6, t7, t12, t13, -t37 * mrSges(7,2) - t56 * t41 + t102; -m(1) * g(2) + t79 * t1 + t76 * t2, t1, t5, t10, t16, t14, -t38 * mrSges(7,3) - t57 * t47 + t101; -m(1) * g(3) + t98, t98, t9, t100, t119, t118, -t54 * mrSges(7,1) + t38 * mrSges(7,2) + t57 * t41 - t69 * t48 + t117;];
f_new  = t3;
