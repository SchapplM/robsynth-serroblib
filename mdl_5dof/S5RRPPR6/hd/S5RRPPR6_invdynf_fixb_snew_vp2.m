% Calculate vector of cutting forces with Newton-Euler
% S5RRPPR6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-31 19:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRPPR6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR6_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR6_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR6_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR6_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR6_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:31:40
% EndTime: 2019-12-31 19:31:45
% DurationCPUTime: 1.99s
% Computational Cost: add. (21242->164), mult. (50063->224), div. (0->0), fcn. (34156->10), ass. (0->84)
t116 = -2 * qJD(3);
t110 = cos(pkin(8));
t106 = qJD(1) * qJD(2);
t86 = sin(qJ(1));
t89 = cos(qJ(1));
t100 = -t89 * g(1) - t86 * g(2);
t91 = qJD(1) ^ 2;
t70 = -t91 * pkin(1) + qJDD(1) * pkin(6) + t100;
t85 = sin(qJ(2));
t111 = t85 * t70;
t112 = pkin(2) * t91;
t88 = cos(qJ(2));
t73 = t85 * qJDD(1) + t88 * t106;
t37 = qJDD(2) * pkin(2) - t73 * qJ(3) - t111 + (qJ(3) * t106 + t85 * t112 - g(3)) * t88;
t102 = -t85 * g(3) + t88 * t70;
t74 = t88 * qJDD(1) - t85 * t106;
t109 = qJD(1) * t85;
t75 = qJD(2) * pkin(2) - qJ(3) * t109;
t80 = t88 ^ 2;
t39 = t74 * qJ(3) - qJD(2) * t75 - t80 * t112 + t102;
t82 = sin(pkin(8));
t67 = (t110 * t85 + t82 * t88) * qJD(1);
t115 = t110 * t37 + t67 * t116 - t82 * t39;
t76 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t109;
t108 = qJD(1) * t88;
t77 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t108;
t66 = -t110 * t108 + t82 * t109;
t104 = t110 * t39 + t66 * t116 + t82 * t37;
t49 = t66 * pkin(3) - t67 * qJ(4);
t90 = qJD(2) ^ 2;
t21 = -t90 * pkin(3) + qJDD(2) * qJ(4) - t66 * t49 + t104;
t53 = -t110 * t74 + t82 * t73;
t54 = t110 * t73 + t82 * t74;
t103 = t86 * g(1) - t89 * g(2);
t97 = -qJDD(1) * pkin(1) - t103;
t94 = -t74 * pkin(2) + qJDD(3) + t75 * t109 + (-qJ(3) * t80 - pkin(6)) * t91 + t97;
t24 = (qJD(2) * t66 - t54) * qJ(4) + (qJD(2) * t67 + t53) * pkin(3) + t94;
t81 = sin(pkin(9));
t83 = cos(pkin(9));
t59 = t81 * qJD(2) + t83 * t67;
t101 = -0.2e1 * qJD(4) * t59 - t81 * t21 + t83 * t24;
t47 = t81 * qJDD(2) + t83 * t54;
t58 = t83 * qJD(2) - t81 * t67;
t15 = (t58 * t66 - t47) * pkin(7) + (t58 * t59 + t53) * pkin(4) + t101;
t105 = 0.2e1 * qJD(4) * t58 + t83 * t21 + t81 * t24;
t45 = t66 * pkin(4) - t59 * pkin(7);
t46 = t83 * qJDD(2) - t81 * t54;
t57 = t58 ^ 2;
t16 = -t57 * pkin(4) + t46 * pkin(7) - t66 * t45 + t105;
t84 = sin(qJ(5));
t87 = cos(qJ(5));
t33 = t87 * t58 - t84 * t59;
t27 = t33 * qJD(5) + t84 * t46 + t87 * t47;
t34 = t84 * t58 + t87 * t59;
t29 = -t33 * mrSges(6,1) + t34 * mrSges(6,2);
t65 = qJD(5) + t66;
t30 = -t65 * mrSges(6,2) + t33 * mrSges(6,3);
t52 = qJDD(5) + t53;
t13 = m(6) * (t87 * t15 - t84 * t16) - t27 * mrSges(6,3) + t52 * mrSges(6,1) - t34 * t29 + t65 * t30;
t26 = -t34 * qJD(5) + t87 * t46 - t84 * t47;
t31 = t65 * mrSges(6,1) - t34 * mrSges(6,3);
t14 = m(6) * (t84 * t15 + t87 * t16) + t26 * mrSges(6,3) - t52 * mrSges(6,2) + t33 * t29 - t65 * t31;
t38 = -t58 * mrSges(5,1) + t59 * mrSges(5,2);
t43 = -t66 * mrSges(5,2) + t58 * mrSges(5,3);
t10 = m(5) * t101 + t53 * mrSges(5,1) - t47 * mrSges(5,3) + t87 * t13 + t84 * t14 - t59 * t38 + t66 * t43;
t44 = t66 * mrSges(5,1) - t59 * mrSges(5,3);
t11 = m(5) * t105 - t53 * mrSges(5,2) + t46 * mrSges(5,3) - t84 * t13 + t87 * t14 + t58 * t38 - t66 * t44;
t60 = -qJD(2) * mrSges(4,2) - t66 * mrSges(4,3);
t61 = qJD(2) * mrSges(4,1) - t67 * mrSges(4,3);
t95 = m(4) * t94 + t53 * mrSges(4,1) + t54 * mrSges(4,2) + t83 * t10 + t81 * t11 + t66 * t60 + t67 * t61;
t114 = (t85 * t76 - t88 * t77) * qJD(1) + m(3) * (-t91 * pkin(6) + t97) - t74 * mrSges(3,1) + t73 * mrSges(3,2) + t95;
t50 = t66 * mrSges(4,1) + t67 * mrSges(4,2);
t20 = -qJDD(2) * pkin(3) - t90 * qJ(4) + t67 * t49 + qJDD(4) - t115;
t96 = t26 * mrSges(6,1) + t33 * t30 - m(6) * (-t46 * pkin(4) - t57 * pkin(7) + t59 * t45 + t20) - t27 * mrSges(6,2) - t34 * t31;
t92 = m(5) * t20 - t46 * mrSges(5,1) + t47 * mrSges(5,2) - t58 * t43 + t59 * t44 - t96;
t12 = m(4) * t115 + qJDD(2) * mrSges(4,1) - t54 * mrSges(4,3) + qJD(2) * t60 - t67 * t50 - t92;
t7 = m(4) * t104 - qJDD(2) * mrSges(4,2) - t53 * mrSges(4,3) - qJD(2) * t61 - t81 * t10 + t83 * t11 - t66 * t50;
t72 = (-mrSges(3,1) * t88 + mrSges(3,2) * t85) * qJD(1);
t4 = m(3) * (-t88 * g(3) - t111) - t73 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t72 * t109 + qJD(2) * t77 + t82 * t7 + t110 * t12;
t5 = m(3) * t102 - qJDD(2) * mrSges(3,2) + t74 * mrSges(3,3) - qJD(2) * t76 + t72 * t108 + t110 * t7 - t82 * t12;
t113 = t88 * t4 + t85 * t5;
t6 = m(2) * t103 + qJDD(1) * mrSges(2,1) - t91 * mrSges(2,2) - t114;
t1 = m(2) * t100 - t91 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t85 * t4 + t88 * t5;
t2 = [-m(1) * g(1) + t89 * t1 - t86 * t6, t1, t5, t7, t11, t14; -m(1) * g(2) + t86 * t1 + t89 * t6, t6, t4, t12, t10, t13; (-m(1) - m(2)) * g(3) + t113, -m(2) * g(3) + t113, t114, t95, t92, -t96;];
f_new = t2;
