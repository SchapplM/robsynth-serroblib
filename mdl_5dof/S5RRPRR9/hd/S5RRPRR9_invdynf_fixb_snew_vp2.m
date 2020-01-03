% Calculate vector of cutting forces with Newton-Euler
% S5RRPRR9
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRPRR9_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR9_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR9_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR9_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR9_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR9_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR9_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:20:07
% EndTime: 2019-12-31 20:20:12
% DurationCPUTime: 1.98s
% Computational Cost: add. (22723->166), mult. (51682->223), div. (0->0), fcn. (35710->10), ass. (0->85)
t106 = qJD(1) * qJD(2);
t86 = sin(qJ(1));
t90 = cos(qJ(1));
t100 = -t90 * g(1) - t86 * g(2);
t92 = qJD(1) ^ 2;
t70 = -t92 * pkin(1) + qJDD(1) * pkin(6) + t100;
t85 = sin(qJ(2));
t111 = t85 * t70;
t112 = pkin(2) * t92;
t89 = cos(qJ(2));
t73 = t85 * qJDD(1) + t89 * t106;
t39 = qJDD(2) * pkin(2) - t73 * qJ(3) - t111 + (qJ(3) * t106 + t85 * t112 - g(3)) * t89;
t103 = -t85 * g(3) + t89 * t70;
t74 = t89 * qJDD(1) - t85 * t106;
t109 = qJD(1) * t85;
t75 = qJD(2) * pkin(2) - qJ(3) * t109;
t80 = t89 ^ 2;
t40 = t74 * qJ(3) - qJD(2) * t75 - t80 * t112 + t103;
t81 = sin(pkin(9));
t82 = cos(pkin(9));
t67 = (t81 * t89 + t82 * t85) * qJD(1);
t115 = -0.2e1 * qJD(3) * t67 + t82 * t39 - t81 * t40;
t76 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t109;
t108 = qJD(1) * t89;
t77 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t108;
t66 = t82 * t108 - t81 * t109;
t105 = 0.2e1 * qJD(3) * t66 + t81 * t39 + t82 * t40;
t50 = -t66 * pkin(3) - t67 * pkin(7);
t91 = qJD(2) ^ 2;
t24 = -t91 * pkin(3) + qJDD(2) * pkin(7) + t66 * t50 + t105;
t54 = -t81 * t73 + t82 * t74;
t55 = t82 * t73 + t81 * t74;
t104 = t86 * g(1) - t90 * g(2);
t98 = -qJDD(1) * pkin(1) - t104;
t95 = -t74 * pkin(2) + qJDD(3) + t75 * t109 + (-qJ(3) * t80 - pkin(6)) * t92 + t98;
t27 = (-qJD(2) * t66 - t55) * pkin(7) + (qJD(2) * t67 - t54) * pkin(3) + t95;
t84 = sin(qJ(4));
t88 = cos(qJ(4));
t102 = -t84 * t24 + t88 * t27;
t57 = t88 * qJD(2) - t84 * t67;
t33 = t57 * qJD(4) + t84 * qJDD(2) + t88 * t55;
t53 = qJDD(4) - t54;
t58 = t84 * qJD(2) + t88 * t67;
t65 = qJD(4) - t66;
t15 = (t57 * t65 - t33) * pkin(8) + (t57 * t58 + t53) * pkin(4) + t102;
t110 = t88 * t24 + t84 * t27;
t32 = -t58 * qJD(4) + t88 * qJDD(2) - t84 * t55;
t47 = t65 * pkin(4) - t58 * pkin(8);
t56 = t57 ^ 2;
t16 = -t56 * pkin(4) + t32 * pkin(8) - t65 * t47 + t110;
t83 = sin(qJ(5));
t87 = cos(qJ(5));
t37 = t87 * t57 - t83 * t58;
t21 = t37 * qJD(5) + t83 * t32 + t87 * t33;
t38 = t83 * t57 + t87 * t58;
t29 = -t37 * mrSges(6,1) + t38 * mrSges(6,2);
t61 = qJD(5) + t65;
t30 = -t61 * mrSges(6,2) + t37 * mrSges(6,3);
t51 = qJDD(5) + t53;
t13 = m(6) * (t87 * t15 - t83 * t16) - t21 * mrSges(6,3) + t51 * mrSges(6,1) - t38 * t29 + t61 * t30;
t20 = -t38 * qJD(5) + t87 * t32 - t83 * t33;
t31 = t61 * mrSges(6,1) - t38 * mrSges(6,3);
t14 = m(6) * (t83 * t15 + t87 * t16) + t20 * mrSges(6,3) - t51 * mrSges(6,2) + t37 * t29 - t61 * t31;
t41 = -t57 * mrSges(5,1) + t58 * mrSges(5,2);
t45 = -t65 * mrSges(5,2) + t57 * mrSges(5,3);
t10 = m(5) * t102 + t53 * mrSges(5,1) - t33 * mrSges(5,3) + t87 * t13 + t83 * t14 - t58 * t41 + t65 * t45;
t46 = t65 * mrSges(5,1) - t58 * mrSges(5,3);
t11 = m(5) * t110 - t53 * mrSges(5,2) + t32 * mrSges(5,3) - t83 * t13 + t87 * t14 + t57 * t41 - t65 * t46;
t59 = -qJD(2) * mrSges(4,2) + t66 * mrSges(4,3);
t60 = qJD(2) * mrSges(4,1) - t67 * mrSges(4,3);
t96 = -m(4) * t95 + t54 * mrSges(4,1) - t55 * mrSges(4,2) - t88 * t10 - t84 * t11 + t66 * t59 - t67 * t60;
t114 = (t85 * t76 - t89 * t77) * qJD(1) + m(3) * (-t92 * pkin(6) + t98) - t74 * mrSges(3,1) + t73 * mrSges(3,2) - t96;
t49 = -t66 * mrSges(4,1) + t67 * mrSges(4,2);
t23 = -qJDD(2) * pkin(3) - t91 * pkin(7) + t67 * t50 - t115;
t97 = t20 * mrSges(6,1) + t37 * t30 - m(6) * (-t32 * pkin(4) - t56 * pkin(8) + t58 * t47 + t23) - t21 * mrSges(6,2) - t38 * t31;
t93 = m(5) * t23 - t32 * mrSges(5,1) + t33 * mrSges(5,2) - t57 * t45 + t58 * t46 - t97;
t12 = m(4) * t115 + qJDD(2) * mrSges(4,1) - t55 * mrSges(4,3) + qJD(2) * t59 - t67 * t49 - t93;
t7 = m(4) * t105 - qJDD(2) * mrSges(4,2) + t54 * mrSges(4,3) - qJD(2) * t60 - t84 * t10 + t88 * t11 + t66 * t49;
t72 = (-mrSges(3,1) * t89 + mrSges(3,2) * t85) * qJD(1);
t4 = m(3) * (-t89 * g(3) - t111) - t73 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t72 * t109 + qJD(2) * t77 + t81 * t7 + t82 * t12;
t5 = m(3) * t103 - qJDD(2) * mrSges(3,2) + t74 * mrSges(3,3) - qJD(2) * t76 + t72 * t108 - t81 * t12 + t82 * t7;
t113 = t89 * t4 + t85 * t5;
t6 = m(2) * t104 + qJDD(1) * mrSges(2,1) - t92 * mrSges(2,2) - t114;
t1 = m(2) * t100 - t92 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t85 * t4 + t89 * t5;
t2 = [-m(1) * g(1) + t90 * t1 - t86 * t6, t1, t5, t7, t11, t14; -m(1) * g(2) + t86 * t1 + t90 * t6, t6, t4, t12, t10, t13; (-m(1) - m(2)) * g(3) + t113, -m(2) * g(3) + t113, t114, -t96, t93, -t97;];
f_new = t2;
