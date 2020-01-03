% Calculate vector of cutting forces with Newton-Euler
% S5RRRRP10
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRRP10_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP10_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP10_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP10_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP10_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP10_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP10_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:08:15
% EndTime: 2019-12-31 22:08:20
% DurationCPUTime: 1.91s
% Computational Cost: add. (22810->170), mult. (48938->228), div. (0->0), fcn. (37244->10), ass. (0->87)
t102 = qJD(1) * qJD(2);
t78 = sin(pkin(5));
t82 = sin(qJ(2));
t86 = cos(qJ(2));
t70 = (-qJDD(1) * t86 + t82 * t102) * t78;
t79 = cos(pkin(5));
t75 = t79 * qJD(1) + qJD(2);
t81 = sin(qJ(3));
t85 = cos(qJ(3));
t105 = qJD(1) * t78;
t98 = t82 * t105;
t58 = t85 * t75 - t81 * t98;
t59 = t81 * t75 + t85 * t98;
t49 = -t58 * pkin(3) - t59 * pkin(9);
t62 = qJDD(3) + t70;
t104 = qJD(1) * t86;
t97 = t78 * t104;
t72 = qJD(3) - t97;
t71 = t72 ^ 2;
t115 = pkin(7) * t78;
t88 = qJD(1) ^ 2;
t83 = sin(qJ(1));
t87 = cos(qJ(1));
t96 = t83 * g(1) - t87 * g(2);
t65 = qJDD(1) * pkin(1) + t88 * t115 + t96;
t113 = t65 * t79;
t93 = -t87 * g(1) - t83 * g(2);
t66 = -t88 * pkin(1) + qJDD(1) * t115 + t93;
t106 = t82 * t113 + t86 * t66;
t68 = (-pkin(2) * t86 - pkin(8) * t82) * t105;
t73 = t75 ^ 2;
t74 = t79 * qJDD(1) + qJDD(2);
t33 = -t73 * pkin(2) + t74 * pkin(8) + (-g(3) * t82 + t68 * t104) * t78 + t106;
t114 = t79 * g(3);
t69 = (qJDD(1) * t82 + t86 * t102) * t78;
t34 = t70 * pkin(2) - t69 * pkin(8) - t114 + (-t65 + (pkin(2) * t82 - pkin(8) * t86) * t75 * qJD(1)) * t78;
t94 = -t81 * t33 + t85 * t34;
t18 = -t62 * pkin(3) - t71 * pkin(9) + t59 * t49 - t94;
t47 = t58 * qJD(3) + t85 * t69 + t81 * t74;
t80 = sin(qJ(4));
t84 = cos(qJ(4));
t52 = t84 * t59 + t80 * t72;
t25 = -t52 * qJD(4) - t80 * t47 + t84 * t62;
t51 = -t80 * t59 + t84 * t72;
t26 = t51 * qJD(4) + t84 * t47 + t80 * t62;
t57 = qJD(4) - t58;
t39 = -t57 * mrSges(6,2) + t51 * mrSges(6,3);
t40 = -t57 * mrSges(5,2) + t51 * mrSges(5,3);
t43 = t57 * mrSges(5,1) - t52 * mrSges(5,3);
t41 = t57 * pkin(4) - t52 * qJ(5);
t42 = t57 * mrSges(6,1) - t52 * mrSges(6,3);
t50 = t51 ^ 2;
t99 = m(6) * (-t25 * pkin(4) - t50 * qJ(5) + t52 * t41 + qJDD(5) + t18) + t26 * mrSges(6,2) + t52 * t42;
t116 = m(5) * t18 + t26 * mrSges(5,2) - (t40 + t39) * t51 - (mrSges(5,1) + mrSges(6,1)) * t25 + t52 * t43 + t99;
t112 = t78 * t82;
t111 = t78 * t86;
t108 = t85 * t33 + t81 * t34;
t19 = -t71 * pkin(3) + t62 * pkin(9) + t58 * t49 + t108;
t92 = -g(3) * t111 + t86 * t113 - t82 * t66;
t32 = -t74 * pkin(2) - t73 * pkin(8) + t68 * t98 - t92;
t46 = -t59 * qJD(3) - t81 * t69 + t85 * t74;
t22 = (-t58 * t72 - t47) * pkin(9) + (t59 * t72 - t46) * pkin(3) + t32;
t109 = t84 * t19 + t80 * t22;
t48 = -t58 * mrSges(4,1) + t59 * mrSges(4,2);
t53 = -t72 * mrSges(4,2) + t58 * mrSges(4,3);
t12 = m(4) * t94 + t62 * mrSges(4,1) - t47 * mrSges(4,3) - t59 * t48 + t72 * t53 - t116;
t63 = t75 * mrSges(3,1) - mrSges(3,3) * t98;
t67 = (-mrSges(3,1) * t86 + mrSges(3,2) * t82) * t105;
t45 = qJDD(4) - t46;
t95 = -t80 * t19 + t84 * t22;
t101 = m(6) * (-0.2e1 * qJD(5) * t52 + (t51 * t57 - t26) * qJ(5) + (t51 * t52 + t45) * pkin(4) + t95) + t57 * t39 + t45 * mrSges(6,1);
t36 = -t51 * mrSges(6,1) + t52 * mrSges(6,2);
t37 = -t51 * mrSges(5,1) + t52 * mrSges(5,2);
t10 = m(5) * t95 + t45 * mrSges(5,1) + t57 * t40 + (-t37 - t36) * t52 + (-mrSges(5,3) - mrSges(6,3)) * t26 + t101;
t100 = m(6) * (-t50 * pkin(4) + t25 * qJ(5) + 0.2e1 * qJD(5) * t51 - t57 * t41 + t109) + t25 * mrSges(6,3) + t51 * t36;
t11 = m(5) * t109 + t25 * mrSges(5,3) + t51 * t37 + (-t43 - t42) * t57 + (-mrSges(5,2) - mrSges(6,2)) * t45 + t100;
t54 = t72 * mrSges(4,1) - t59 * mrSges(4,3);
t9 = m(4) * t108 - t62 * mrSges(4,2) + t46 * mrSges(4,3) - t80 * t10 + t84 * t11 + t58 * t48 - t72 * t54;
t4 = m(3) * (-g(3) * t112 + t106) - t70 * mrSges(3,3) - t74 * mrSges(3,2) + t67 * t97 - t75 * t63 + t85 * t9 - t81 * t12;
t64 = -t75 * mrSges(3,2) + mrSges(3,3) * t97;
t6 = m(3) * (-t78 * t65 - t114) + t69 * mrSges(3,2) + t70 * mrSges(3,1) + t81 * t9 + t85 * t12 + (t63 * t82 - t64 * t86) * t105;
t89 = m(4) * t32 - t46 * mrSges(4,1) + t47 * mrSges(4,2) + t84 * t10 + t80 * t11 - t58 * t53 + t59 * t54;
t8 = m(3) * t92 + t74 * mrSges(3,1) - t69 * mrSges(3,3) + t75 * t64 - t67 * t98 - t89;
t103 = t8 * t111 + t4 * t112 + t79 * t6;
t2 = m(2) * t93 - t88 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t86 * t4 - t82 * t8;
t1 = m(2) * t96 + qJDD(1) * mrSges(2,1) - t88 * mrSges(2,2) - t78 * t6 + (t82 * t4 + t86 * t8) * t79;
t3 = [-m(1) * g(1) - t83 * t1 + t87 * t2, t2, t4, t9, t11, -t45 * mrSges(6,2) - t57 * t42 + t100; -m(1) * g(2) + t87 * t1 + t83 * t2, t1, t8, t12, t10, -t26 * mrSges(6,3) - t52 * t36 + t101; (-m(1) - m(2)) * g(3) + t103, -m(2) * g(3) + t103, t6, t89, t116, -t25 * mrSges(6,1) - t51 * t39 + t99;];
f_new = t3;
