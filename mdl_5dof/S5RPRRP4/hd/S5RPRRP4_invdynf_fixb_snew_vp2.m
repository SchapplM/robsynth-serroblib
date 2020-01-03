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
% Datum: 2020-01-03 11:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:49:19
% EndTime: 2020-01-03 11:49:22
% DurationCPUTime: 1.10s
% Computational Cost: add. (8562->150), mult. (20533->196), div. (0->0), fcn. (13374->8), ass. (0->79)
t75 = sin(pkin(8));
t105 = qJD(1) * t75;
t77 = sin(qJ(4));
t78 = sin(qJ(3));
t80 = cos(qJ(4));
t81 = cos(qJ(3));
t47 = (-t77 * t78 + t80 * t81) * t105;
t102 = qJD(1) * qJD(3);
t54 = (-qJDD(1) * t78 - t81 * t102) * t75;
t55 = (qJDD(1) * t81 - t78 * t102) * t75;
t26 = -t47 * qJD(4) + t80 * t54 - t77 * t55;
t46 = (-t77 * t81 - t78 * t80) * t105;
t27 = t46 * qJD(4) + t77 * t54 + t80 * t55;
t76 = cos(pkin(8));
t104 = t76 * qJD(1);
t65 = qJD(3) - t104;
t63 = qJD(4) + t65;
t36 = -t63 * mrSges(6,2) + t46 * mrSges(6,3);
t37 = -t63 * mrSges(5,2) + t46 * mrSges(5,3);
t40 = t63 * mrSges(5,1) - t47 * mrSges(5,3);
t72 = t75 ^ 2;
t83 = qJD(1) ^ 2;
t113 = t72 * t83;
t101 = t78 ^ 2 * t113;
t79 = sin(qJ(1));
t82 = cos(qJ(1));
t94 = -t79 * g(2) + t82 * g(3);
t58 = -t83 * pkin(1) + qJDD(1) * qJ(2) + t94;
t108 = -t76 * g(1) - t75 * t58;
t91 = -pkin(2) * t76 - pkin(6) * t75;
t60 = t91 * qJD(1);
t95 = 0.2e1 * qJD(1) * qJD(2);
t31 = t60 * t105 + t75 * t95 - t108;
t96 = t81 * t105;
t53 = t65 * pkin(3) - pkin(7) * t96;
t85 = -t54 * pkin(3) - pkin(7) * t101 + t53 * t96 + t31;
t38 = t63 * pkin(4) - t47 * qJ(5);
t39 = t63 * mrSges(6,1) - t47 * mrSges(6,3);
t45 = t46 ^ 2;
t98 = m(6) * (-t26 * pkin(4) - t45 * qJ(5) + t47 * t38 + qJDD(5) + t85) + t27 * mrSges(6,2) + t47 * t39;
t118 = -m(5) * t85 - t27 * mrSges(5,2) + (t37 + t36) * t46 + (mrSges(5,1) + mrSges(6,1)) * t26 - t47 * t40 - t98;
t121 = -m(4) * t31 + t54 * mrSges(4,1) - t55 * mrSges(4,2) + t118;
t120 = (t76 ^ 2 + t72) * mrSges(3,3);
t103 = qJDD(1) * mrSges(3,3);
t90 = -t76 * mrSges(3,1) + t75 * mrSges(3,2);
t59 = t90 * qJD(1);
t97 = t78 * t105;
t50 = -t65 * mrSges(4,2) - mrSges(4,3) * t97;
t51 = t65 * mrSges(4,1) - mrSges(4,3) * t96;
t89 = t50 * t78 + t51 * t81;
t10 = m(3) * t108 + (-t103 + (-0.2e1 * m(3) * qJD(2) - t59 - t89) * qJD(1)) * t75 + t121;
t92 = -t75 * g(1) + (t58 + t95) * t76;
t32 = t60 * t104 + t92;
t107 = -t82 * g(2) - t79 * g(3);
t87 = -t83 * qJ(2) + qJDD(2) - t107;
t43 = (-pkin(1) + t91) * qJDD(1) + t87;
t42 = t81 * t43;
t52 = (mrSges(4,1) * t78 + mrSges(4,2) * t81) * t105;
t64 = -t76 * qJDD(1) + qJDD(3);
t62 = qJDD(4) + t64;
t18 = t64 * pkin(3) - t55 * pkin(7) + t42 + (-pkin(3) * t81 * t113 - pkin(7) * t65 * t105 - t32) * t78;
t110 = t81 * t32 + t78 * t43;
t19 = -pkin(3) * t101 + t54 * pkin(7) - t65 * t53 + t110;
t93 = t80 * t18 - t77 * t19;
t100 = m(6) * (-0.2e1 * qJD(5) * t47 + (t46 * t63 - t27) * qJ(5) + (t46 * t47 + t62) * pkin(4) + t93) + t63 * t36 + t62 * mrSges(6,1);
t33 = -t46 * mrSges(6,1) + t47 * mrSges(6,2);
t34 = -t46 * mrSges(5,1) + t47 * mrSges(5,2);
t7 = m(5) * t93 + t62 * mrSges(5,1) + t63 * t37 + (-t34 - t33) * t47 + (-mrSges(5,3) - mrSges(6,3)) * t27 + t100;
t111 = t77 * t18 + t80 * t19;
t99 = m(6) * (-t45 * pkin(4) + t26 * qJ(5) + 0.2e1 * qJD(5) * t46 - t63 * t38 + t111) + t46 * t33 + t26 * mrSges(6,3);
t8 = m(5) * t111 + t26 * mrSges(5,3) + t46 * t34 + (-t40 - t39) * t63 + (-mrSges(5,2) - mrSges(6,2)) * t62 + t99;
t5 = m(4) * (-t78 * t32 + t42) - t55 * mrSges(4,3) + t64 * mrSges(4,1) - t52 * t96 + t65 * t50 + t77 * t8 + t80 * t7;
t6 = m(4) * t110 - t64 * mrSges(4,2) + t54 * mrSges(4,3) - t65 * t51 - t52 * t97 - t77 * t7 + t80 * t8;
t4 = m(3) * t92 + t81 * t6 - t78 * t5 + (qJD(1) * t59 + t103) * t76;
t117 = t76 * t10 + t75 * t4;
t86 = m(3) * (-qJDD(1) * pkin(1) + t87) + t81 * t5 + t78 * t6;
t2 = m(2) * t107 + (-mrSges(2,2) + t120) * t83 + (mrSges(2,1) - t90) * qJDD(1) - t86;
t1 = m(2) * t94 - t83 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t75 * t10 + t76 * t4;
t3 = [(-m(1) - m(2)) * g(1) + t117, t1, t4, t6, t8, -t62 * mrSges(6,2) - t63 * t39 + t99; -m(1) * g(2) + t79 * t1 + t82 * t2, t2, t10, t5, t7, -t27 * mrSges(6,3) - t47 * t33 + t100; -m(1) * g(3) - t82 * t1 + t79 * t2, -m(2) * g(1) + t117, t90 * qJDD(1) - t83 * t120 + t86, t89 * t105 - t121, -t118, -t26 * mrSges(6,1) - t46 * t36 + t98;];
f_new = t3;
