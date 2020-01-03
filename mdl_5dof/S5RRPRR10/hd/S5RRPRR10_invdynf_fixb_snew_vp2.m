% Calculate vector of cutting forces with Newton-Euler
% S5RRPRR10
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRPRR10_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR10_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR10_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR10_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR10_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR10_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR10_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:23:58
% EndTime: 2019-12-31 20:24:06
% DurationCPUTime: 3.24s
% Computational Cost: add. (36805->171), mult. (96541->244), div. (0->0), fcn. (74059->12), ass. (0->95)
t122 = -2 * qJD(3);
t83 = sin(pkin(5));
t114 = qJD(1) * t83;
t82 = sin(pkin(10));
t84 = cos(pkin(10));
t88 = sin(qJ(2));
t92 = cos(qJ(2));
t65 = (t82 * t88 - t84 * t92) * t114;
t89 = sin(qJ(1));
t93 = cos(qJ(1));
t106 = t89 * g(1) - t93 * g(2);
t120 = pkin(7) * t83;
t94 = qJD(1) ^ 2;
t71 = qJDD(1) * pkin(1) + t94 * t120 + t106;
t85 = cos(pkin(5));
t119 = t71 * t85;
t103 = -t93 * g(1) - t89 * g(2);
t72 = -t94 * pkin(1) + qJDD(1) * t120 + t103;
t104 = t92 * t119 - t88 * t72;
t118 = t83 ^ 2 * t94;
t111 = qJD(1) * qJD(2);
t74 = (qJDD(1) * t88 + t92 * t111) * t83;
t77 = t85 * qJDD(1) + qJDD(2);
t78 = t85 * qJD(1) + qJD(2);
t31 = t77 * pkin(2) - t74 * qJ(3) + (pkin(2) * t88 * t118 + (qJ(3) * qJD(1) * t78 - g(3)) * t83) * t92 + t104;
t110 = t92 ^ 2 * t118;
t108 = t88 * t114;
t68 = t78 * pkin(2) - qJ(3) * t108;
t75 = (qJDD(1) * t92 - t88 * t111) * t83;
t117 = t83 * t88;
t99 = -g(3) * t117 + t88 * t119 + t92 * t72;
t34 = -pkin(2) * t110 + t75 * qJ(3) - t78 * t68 + t99;
t66 = (t82 * t92 + t84 * t88) * t114;
t121 = t66 * t122 + t84 * t31 - t82 * t34;
t116 = t83 * t92;
t109 = t65 * t122 + t82 * t31 + t84 * t34;
t48 = t65 * pkin(3) - t66 * pkin(8);
t76 = t78 ^ 2;
t22 = -t76 * pkin(3) + t77 * pkin(8) - t65 * t48 + t109;
t51 = -t82 * t74 + t84 * t75;
t52 = t84 * t74 + t82 * t75;
t102 = -t85 * g(3) - t83 * t71;
t96 = -t75 * pkin(2) - qJ(3) * t110 + t68 * t108 + qJDD(3) + t102;
t24 = (t65 * t78 - t52) * pkin(8) + (t66 * t78 - t51) * pkin(3) + t96;
t87 = sin(qJ(4));
t91 = cos(qJ(4));
t115 = t91 * t22 + t87 * t24;
t47 = t65 * mrSges(4,1) + t66 * mrSges(4,2);
t56 = -t78 * mrSges(4,2) - t65 * mrSges(4,3);
t54 = -t87 * t66 + t91 * t78;
t55 = t91 * t66 + t87 * t78;
t40 = -t54 * pkin(4) - t55 * pkin(9);
t50 = qJDD(4) - t51;
t64 = qJD(4) + t65;
t63 = t64 ^ 2;
t18 = -t63 * pkin(4) + t50 * pkin(9) + t54 * t40 + t115;
t21 = -t77 * pkin(3) - t76 * pkin(8) + t66 * t48 - t121;
t36 = -t55 * qJD(4) - t87 * t52 + t91 * t77;
t37 = t54 * qJD(4) + t91 * t52 + t87 * t77;
t19 = (-t54 * t64 - t37) * pkin(9) + (t55 * t64 - t36) * pkin(4) + t21;
t86 = sin(qJ(5));
t90 = cos(qJ(5));
t42 = -t86 * t55 + t90 * t64;
t26 = t42 * qJD(5) + t90 * t37 + t86 * t50;
t43 = t90 * t55 + t86 * t64;
t27 = -t42 * mrSges(6,1) + t43 * mrSges(6,2);
t53 = qJD(5) - t54;
t32 = -t53 * mrSges(6,2) + t42 * mrSges(6,3);
t35 = qJDD(5) - t36;
t15 = m(6) * (-t86 * t18 + t90 * t19) - t26 * mrSges(6,3) + t35 * mrSges(6,1) - t43 * t27 + t53 * t32;
t25 = -t43 * qJD(5) - t86 * t37 + t90 * t50;
t33 = t53 * mrSges(6,1) - t43 * mrSges(6,3);
t16 = m(6) * (t90 * t18 + t86 * t19) + t25 * mrSges(6,3) - t35 * mrSges(6,2) + t42 * t27 - t53 * t33;
t44 = -t64 * mrSges(5,2) + t54 * mrSges(5,3);
t45 = t64 * mrSges(5,1) - t55 * mrSges(5,3);
t95 = m(5) * t21 - t36 * mrSges(5,1) + t37 * mrSges(5,2) + t90 * t15 + t86 * t16 - t54 * t44 + t55 * t45;
t10 = m(4) * t121 + t77 * mrSges(4,1) - t52 * mrSges(4,3) - t66 * t47 + t78 * t56 - t95;
t107 = (-mrSges(3,1) * t92 + mrSges(3,2) * t88) * t114 ^ 2;
t39 = -t54 * mrSges(5,1) + t55 * mrSges(5,2);
t12 = m(5) * t115 - t50 * mrSges(5,2) + t36 * mrSges(5,3) - t86 * t15 + t90 * t16 + t54 * t39 - t64 * t45;
t101 = -t87 * t22 + t91 * t24;
t97 = m(6) * (-t50 * pkin(4) - t63 * pkin(9) + t55 * t40 - t101) - t25 * mrSges(6,1) + t26 * mrSges(6,2) - t42 * t32 + t43 * t33;
t14 = m(5) * t101 + t50 * mrSges(5,1) - t37 * mrSges(5,3) - t55 * t39 + t64 * t44 - t97;
t57 = t78 * mrSges(4,1) - t66 * mrSges(4,3);
t7 = m(4) * t109 - t77 * mrSges(4,2) + t51 * mrSges(4,3) + t91 * t12 - t87 * t14 - t65 * t47 - t78 * t57;
t70 = t92 * mrSges(3,3) * t114 - t78 * mrSges(3,2);
t5 = m(3) * (-g(3) * t116 + t104) - t74 * mrSges(3,3) + t77 * mrSges(3,1) - t88 * t107 + t78 * t70 + t82 * t7 + t84 * t10;
t69 = t78 * mrSges(3,1) - mrSges(3,3) * t108;
t6 = m(3) * t99 - t77 * mrSges(3,2) + t75 * mrSges(3,3) - t82 * t10 + t92 * t107 - t78 * t69 + t84 * t7;
t98 = m(4) * t96 - t51 * mrSges(4,1) + t52 * mrSges(4,2) + t87 * t12 + t91 * t14 + t65 * t56 + t66 * t57;
t9 = m(3) * t102 + t74 * mrSges(3,2) - t75 * mrSges(3,1) + (t69 * t88 - t70 * t92) * t114 + t98;
t112 = t5 * t116 + t6 * t117 + t85 * t9;
t2 = m(2) * t103 - t94 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t88 * t5 + t92 * t6;
t1 = m(2) * t106 + qJDD(1) * mrSges(2,1) - t94 * mrSges(2,2) - t83 * t9 + (t92 * t5 + t88 * t6) * t85;
t3 = [-m(1) * g(1) - t89 * t1 + t93 * t2, t2, t6, t7, t12, t16; -m(1) * g(2) + t93 * t1 + t89 * t2, t1, t5, t10, t14, t15; (-m(1) - m(2)) * g(3) + t112, -m(2) * g(3) + t112, t9, t98, t95, t97;];
f_new = t3;
