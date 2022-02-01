% Calculate vector of cutting forces with Newton-Euler
% S5RPRRP4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2022-01-23 09:33
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRRP4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP4_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP4_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP4_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:31:54
% EndTime: 2022-01-23 09:31:56
% DurationCPUTime: 1.09s
% Computational Cost: add. (8562->150), mult. (20533->196), div. (0->0), fcn. (13374->8), ass. (0->79)
t73 = sin(pkin(8));
t104 = qJD(1) * t73;
t75 = sin(qJ(4));
t76 = sin(qJ(3));
t78 = cos(qJ(4));
t79 = cos(qJ(3));
t47 = (-t75 * t76 + t78 * t79) * t104;
t101 = qJD(1) * qJD(3);
t54 = (-qJDD(1) * t76 - t79 * t101) * t73;
t55 = (qJDD(1) * t79 - t76 * t101) * t73;
t26 = -t47 * qJD(4) + t78 * t54 - t75 * t55;
t46 = (-t75 * t79 - t76 * t78) * t104;
t27 = t46 * qJD(4) + t75 * t54 + t78 * t55;
t74 = cos(pkin(8));
t103 = t74 * qJD(1);
t65 = qJD(3) - t103;
t63 = qJD(4) + t65;
t36 = -t63 * mrSges(6,2) + t46 * mrSges(6,3);
t37 = -t63 * mrSges(5,2) + t46 * mrSges(5,3);
t40 = t63 * mrSges(5,1) - t47 * mrSges(5,3);
t70 = t73 ^ 2;
t81 = qJD(1) ^ 2;
t111 = t70 * t81;
t100 = t76 ^ 2 * t111;
t77 = sin(qJ(1));
t80 = cos(qJ(1));
t89 = -t80 * g(1) - t77 * g(2);
t58 = -t81 * pkin(1) + qJDD(1) * qJ(2) + t89;
t106 = -t74 * g(3) - t73 * t58;
t90 = -pkin(2) * t74 - pkin(6) * t73;
t60 = t90 * qJD(1);
t94 = 0.2e1 * qJD(1) * qJD(2);
t31 = t60 * t104 + t73 * t94 - t106;
t95 = t79 * t104;
t53 = t65 * pkin(3) - pkin(7) * t95;
t83 = -t54 * pkin(3) - pkin(7) * t100 + t53 * t95 + t31;
t38 = t63 * pkin(4) - t47 * qJ(5);
t39 = t63 * mrSges(6,1) - t47 * mrSges(6,3);
t45 = t46 ^ 2;
t97 = m(6) * (-t26 * pkin(4) - t45 * qJ(5) + t47 * t38 + qJDD(5) + t83) + t27 * mrSges(6,2) + t47 * t39;
t116 = -m(5) * t83 - t27 * mrSges(5,2) + (t37 + t36) * t46 + (mrSges(5,1) + mrSges(6,1)) * t26 - t47 * t40 - t97;
t119 = -m(4) * t31 + t54 * mrSges(4,1) - t55 * mrSges(4,2) + t116;
t118 = (t74 ^ 2 + t70) * mrSges(3,3);
t102 = qJDD(1) * mrSges(3,3);
t88 = -t74 * mrSges(3,1) + t73 * mrSges(3,2);
t59 = t88 * qJD(1);
t96 = t76 * t104;
t50 = -t65 * mrSges(4,2) - mrSges(4,3) * t96;
t51 = t65 * mrSges(4,1) - mrSges(4,3) * t95;
t87 = t50 * t76 + t51 * t79;
t10 = m(3) * t106 + (-t102 + (-0.2e1 * m(3) * qJD(2) - t59 - t87) * qJD(1)) * t73 + t119;
t91 = -t73 * g(3) + (t58 + t94) * t74;
t32 = t60 * t103 + t91;
t93 = t77 * g(1) - t80 * g(2);
t84 = -t81 * qJ(2) + qJDD(2) - t93;
t43 = (-pkin(1) + t90) * qJDD(1) + t84;
t42 = t79 * t43;
t52 = (mrSges(4,1) * t76 + mrSges(4,2) * t79) * t104;
t64 = -t74 * qJDD(1) + qJDD(3);
t33 = -t46 * mrSges(6,1) + t47 * mrSges(6,2);
t34 = -t46 * mrSges(5,1) + t47 * mrSges(5,2);
t62 = qJDD(4) + t64;
t18 = t64 * pkin(3) - t55 * pkin(7) + t42 + (-pkin(3) * t79 * t111 - pkin(7) * t65 * t104 - t32) * t76;
t108 = t79 * t32 + t76 * t43;
t19 = -pkin(3) * t100 + t54 * pkin(7) - t65 * t53 + t108;
t92 = t78 * t18 - t75 * t19;
t99 = m(6) * (-0.2e1 * qJD(5) * t47 + (t46 * t63 - t27) * qJ(5) + (t46 * t47 + t62) * pkin(4) + t92) + t63 * t36 + t62 * mrSges(6,1);
t7 = m(5) * t92 + t62 * mrSges(5,1) + t63 * t37 + (-t34 - t33) * t47 + (-mrSges(5,3) - mrSges(6,3)) * t27 + t99;
t109 = t75 * t18 + t78 * t19;
t98 = m(6) * (-t45 * pkin(4) + t26 * qJ(5) + 0.2e1 * qJD(5) * t46 - t63 * t38 + t109) + t46 * t33 + t26 * mrSges(6,3);
t8 = m(5) * t109 + t26 * mrSges(5,3) + t46 * t34 + (-t40 - t39) * t63 + (-mrSges(5,2) - mrSges(6,2)) * t62 + t98;
t5 = m(4) * (-t76 * t32 + t42) - t55 * mrSges(4,3) + t64 * mrSges(4,1) - t52 * t95 + t65 * t50 + t75 * t8 + t78 * t7;
t6 = m(4) * t108 - t64 * mrSges(4,2) + t54 * mrSges(4,3) - t65 * t51 - t52 * t96 - t75 * t7 + t78 * t8;
t4 = m(3) * t91 + t79 * t6 - t76 * t5 + (qJD(1) * t59 + t102) * t74;
t115 = t74 * t10 + t73 * t4;
t85 = m(3) * (-qJDD(1) * pkin(1) + t84) + t79 * t5 + t76 * t6;
t2 = m(2) * t93 + (-mrSges(2,2) + t118) * t81 + (mrSges(2,1) - t88) * qJDD(1) - t85;
t1 = m(2) * t89 - t81 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t73 * t10 + t74 * t4;
t3 = [-m(1) * g(1) + t80 * t1 - t77 * t2, t1, t4, t6, t8, -t62 * mrSges(6,2) - t63 * t39 + t98; -m(1) * g(2) + t77 * t1 + t80 * t2, t2, t10, t5, t7, -t27 * mrSges(6,3) - t47 * t33 + t99; (-m(1) - m(2)) * g(3) + t115, -m(2) * g(3) + t115, t88 * qJDD(1) - t81 * t118 + t85, t87 * t104 - t119, -t116, -t26 * mrSges(6,1) - t46 * t36 + t97;];
f_new = t3;
