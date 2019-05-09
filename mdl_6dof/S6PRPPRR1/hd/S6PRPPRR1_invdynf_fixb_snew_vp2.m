% Calculate vector of cutting forces with Newton-Euler
% S6PRPPRR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
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
% Datum: 2019-05-04 21:43
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRPPRR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR1_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPPRR1_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR1_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR1_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPPRR1_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 21:39:54
% EndTime: 2019-05-04 21:39:58
% DurationCPUTime: 2.36s
% Computational Cost: add. (32781->139), mult. (66510->188), div. (0->0), fcn. (49112->14), ass. (0->89)
t91 = qJD(2) ^ 2;
t80 = cos(pkin(12));
t74 = t80 ^ 2;
t76 = sin(pkin(12));
t113 = t76 ^ 2 + t74;
t121 = t113 * mrSges(5,3);
t85 = sin(qJ(5));
t88 = cos(qJ(5));
t97 = t76 * t85 - t80 * t88;
t59 = t97 * qJD(2);
t98 = t76 * t88 + t80 * t85;
t60 = t98 * qJD(2);
t110 = t60 * qJD(5);
t48 = -t97 * qJDD(2) - t110;
t120 = pkin(4) * t91;
t78 = sin(pkin(10));
t82 = cos(pkin(10));
t65 = t78 * g(1) - t82 * g(2);
t83 = cos(pkin(6));
t119 = t65 * t83;
t79 = sin(pkin(6));
t86 = sin(qJ(2));
t118 = t79 * t86;
t89 = cos(qJ(2));
t117 = t79 * t89;
t112 = pkin(8) * qJDD(2);
t109 = qJD(2) * qJD(4);
t75 = -g(3) + qJDD(1);
t103 = -t79 * t65 + t83 * t75;
t55 = qJDD(3) + t103;
t114 = -0.2e1 * t76 * t109 + t80 * t55;
t66 = -t82 * g(1) - t78 * g(2);
t102 = t75 * t117 + t89 * t119 - t86 * t66;
t41 = qJDD(2) * pkin(2) + t102;
t106 = t75 * t118 + t86 * t119 + t89 * t66;
t42 = -t91 * pkin(2) + t106;
t77 = sin(pkin(11));
t81 = cos(pkin(11));
t115 = t77 * t41 + t81 * t42;
t32 = -t91 * pkin(3) + qJDD(2) * qJ(4) + t115;
t26 = (t80 * t120 - t112 - t32) * t76 + t114;
t108 = t76 * t55 + (0.2e1 * t109 + t32) * t80;
t27 = t80 * t112 - t74 * t120 + t108;
t116 = t85 * t26 + t88 * t27;
t111 = t59 * qJD(5);
t47 = t59 * pkin(5) - t60 * pkin(9);
t90 = qJD(5) ^ 2;
t22 = -t90 * pkin(5) + qJDD(5) * pkin(9) - t59 * t47 + t116;
t49 = t98 * qJDD(2) - t111;
t104 = t81 * t41 - t77 * t42;
t101 = qJDD(4) - t104;
t92 = (-pkin(4) * t80 - pkin(3)) * qJDD(2) + (-t113 * pkin(8) - qJ(4)) * t91 + t101;
t23 = (-t49 + t111) * pkin(9) + (-t48 + t110) * pkin(5) + t92;
t84 = sin(qJ(6));
t87 = cos(qJ(6));
t52 = t87 * qJD(5) - t84 * t60;
t34 = t52 * qJD(6) + t84 * qJDD(5) + t87 * t49;
t53 = t84 * qJD(5) + t87 * t60;
t35 = -t52 * mrSges(7,1) + t53 * mrSges(7,2);
t58 = qJD(6) + t59;
t36 = -t58 * mrSges(7,2) + t52 * mrSges(7,3);
t46 = qJDD(6) - t48;
t19 = m(7) * (-t84 * t22 + t87 * t23) - t34 * mrSges(7,3) + t46 * mrSges(7,1) - t53 * t35 + t58 * t36;
t33 = -t53 * qJD(6) + t87 * qJDD(5) - t84 * t49;
t37 = t58 * mrSges(7,1) - t53 * mrSges(7,3);
t20 = m(7) * (t87 * t22 + t84 * t23) + t33 * mrSges(7,3) - t46 * mrSges(7,2) + t52 * t35 - t58 * t37;
t44 = t59 * mrSges(6,1) + t60 * mrSges(6,2);
t57 = qJD(5) * mrSges(6,1) - t60 * mrSges(6,3);
t15 = m(6) * t116 - qJDD(5) * mrSges(6,2) + t48 * mrSges(6,3) - qJD(5) * t57 - t84 * t19 + t87 * t20 - t59 * t44;
t56 = -qJD(5) * mrSges(6,2) - t59 * mrSges(6,3);
t99 = t88 * t26 - t85 * t27;
t93 = m(7) * (-qJDD(5) * pkin(5) - t90 * pkin(9) + t60 * t47 - t99) - t33 * mrSges(7,1) + t34 * mrSges(7,2) - t52 * t36 + t53 * t37;
t16 = m(6) * t99 + qJDD(5) * mrSges(6,1) - t49 * mrSges(6,3) + qJD(5) * t56 - t60 * t44 - t93;
t100 = -t80 * mrSges(5,1) + t76 * mrSges(5,2);
t96 = qJDD(2) * mrSges(5,3) + t91 * t100;
t12 = m(5) * t114 + t85 * t15 + t88 * t16 + (-m(5) * t32 - t96) * t76;
t13 = m(5) * t108 + t88 * t15 - t85 * t16 + t96 * t80;
t107 = m(4) * t55 + t80 * t12 + t76 * t13;
t95 = -m(6) * t92 + t48 * mrSges(6,1) - t49 * mrSges(6,2) - t87 * t19 - t84 * t20 - t59 * t56 - t60 * t57;
t94 = m(5) * (-qJDD(2) * pkin(3) - t91 * qJ(4) + t101) - t95;
t14 = m(4) * t104 + (-mrSges(4,2) + t121) * t91 + (mrSges(4,1) - t100) * qJDD(2) - t94;
t7 = m(4) * t115 - t91 * mrSges(4,1) - qJDD(2) * mrSges(4,2) - t76 * t12 + t80 * t13;
t5 = m(3) * t102 + qJDD(2) * mrSges(3,1) - t91 * mrSges(3,2) + t81 * t14 + t77 * t7;
t6 = m(3) * t106 - t91 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t77 * t14 + t81 * t7;
t9 = m(3) * t103 + t107;
t105 = m(2) * t75 + t5 * t117 + t6 * t118 + t83 * t9;
t2 = m(2) * t66 - t86 * t5 + t89 * t6;
t1 = m(2) * t65 - t79 * t9 + (t5 * t89 + t6 * t86) * t83;
t3 = [-m(1) * g(1) - t78 * t1 + t82 * t2, t2, t6, t7, t13, t15, t20; -m(1) * g(2) + t82 * t1 + t78 * t2, t1, t5, t14, t12, t16, t19; -m(1) * g(3) + t105, t105, t9, t107, t100 * qJDD(2) - t91 * t121 + t94, -t95, t93;];
f_new  = t3;
