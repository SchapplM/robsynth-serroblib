% Calculate vector of cutting forces with Newton-Euler
% S5RRRRR8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRRR8_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR8_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR8_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR8_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR8_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR8_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR8_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR8_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR8_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:24:28
% EndTime: 2019-12-31 22:24:33
% DurationCPUTime: 2.16s
% Computational Cost: add. (27599->167), mult. (55543->221), div. (0->0), fcn. (38998->10), ass. (0->87)
t106 = qJD(1) * qJD(2);
t86 = sin(qJ(2));
t91 = cos(qJ(2));
t72 = t86 * qJDD(1) + t91 * t106;
t73 = t91 * qJDD(1) - t86 * t106;
t108 = qJD(1) * t86;
t74 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t108;
t107 = qJD(1) * t91;
t75 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t107;
t93 = qJD(1) ^ 2;
t85 = sin(qJ(3));
t90 = cos(qJ(3));
t66 = (t85 * t91 + t86 * t90) * qJD(1);
t47 = -t66 * qJD(3) - t85 * t72 + t90 * t73;
t65 = t90 * t107 - t85 * t108;
t48 = t65 * qJD(3) + t90 * t72 + t85 * t73;
t81 = qJD(2) + qJD(3);
t76 = qJD(2) * pkin(2) - pkin(7) * t108;
t82 = t91 ^ 2;
t87 = sin(qJ(1));
t92 = cos(qJ(1));
t105 = t87 * g(1) - t92 * g(2);
t99 = -qJDD(1) * pkin(1) - t105;
t96 = -t73 * pkin(2) + t76 * t108 + (-pkin(7) * t82 - pkin(6)) * t93 + t99;
t24 = (-t65 * t81 - t48) * pkin(8) + (t66 * t81 - t47) * pkin(3) + t96;
t101 = -t92 * g(1) - t87 * g(2);
t68 = -t93 * pkin(1) + qJDD(1) * pkin(6) + t101;
t111 = t86 * t68;
t112 = pkin(2) * t93;
t40 = qJDD(2) * pkin(2) - t72 * pkin(7) - t111 + (pkin(7) * t106 + t86 * t112 - g(3)) * t91;
t104 = -t86 * g(3) + t91 * t68;
t41 = t73 * pkin(7) - qJD(2) * t76 - t82 * t112 + t104;
t109 = t85 * t40 + t90 * t41;
t55 = -t65 * pkin(3) - t66 * pkin(8);
t79 = t81 ^ 2;
t80 = qJDD(2) + qJDD(3);
t27 = -t79 * pkin(3) + t80 * pkin(8) + t65 * t55 + t109;
t84 = sin(qJ(4));
t89 = cos(qJ(4));
t103 = t89 * t24 - t84 * t27;
t57 = -t84 * t66 + t89 * t81;
t31 = t57 * qJD(4) + t89 * t48 + t84 * t80;
t45 = qJDD(4) - t47;
t58 = t89 * t66 + t84 * t81;
t64 = qJD(4) - t65;
t15 = (t57 * t64 - t31) * pkin(9) + (t57 * t58 + t45) * pkin(4) + t103;
t110 = t84 * t24 + t89 * t27;
t30 = -t58 * qJD(4) - t84 * t48 + t89 * t80;
t52 = t64 * pkin(4) - t58 * pkin(9);
t56 = t57 ^ 2;
t16 = -t56 * pkin(4) + t30 * pkin(9) - t64 * t52 + t110;
t83 = sin(qJ(5));
t88 = cos(qJ(5));
t34 = t88 * t57 - t83 * t58;
t21 = t34 * qJD(5) + t83 * t30 + t88 * t31;
t35 = t83 * t57 + t88 * t58;
t29 = -t34 * mrSges(6,1) + t35 * mrSges(6,2);
t62 = qJD(5) + t64;
t32 = -t62 * mrSges(6,2) + t34 * mrSges(6,3);
t42 = qJDD(5) + t45;
t13 = m(6) * (t88 * t15 - t83 * t16) - t21 * mrSges(6,3) + t42 * mrSges(6,1) - t35 * t29 + t62 * t32;
t20 = -t35 * qJD(5) + t88 * t30 - t83 * t31;
t33 = t62 * mrSges(6,1) - t35 * mrSges(6,3);
t14 = m(6) * (t83 * t15 + t88 * t16) + t20 * mrSges(6,3) - t42 * mrSges(6,2) + t34 * t29 - t62 * t33;
t39 = -t57 * mrSges(5,1) + t58 * mrSges(5,2);
t50 = -t64 * mrSges(5,2) + t57 * mrSges(5,3);
t10 = m(5) * t103 + t45 * mrSges(5,1) - t31 * mrSges(5,3) + t88 * t13 + t83 * t14 - t58 * t39 + t64 * t50;
t51 = t64 * mrSges(5,1) - t58 * mrSges(5,3);
t11 = m(5) * t110 - t45 * mrSges(5,2) + t30 * mrSges(5,3) - t83 * t13 + t88 * t14 + t57 * t39 - t64 * t51;
t59 = -t81 * mrSges(4,2) + t65 * mrSges(4,3);
t60 = t81 * mrSges(4,1) - t66 * mrSges(4,3);
t97 = -m(4) * t96 + t47 * mrSges(4,1) - t48 * mrSges(4,2) - t89 * t10 - t84 * t11 + t65 * t59 - t66 * t60;
t114 = (t86 * t74 - t91 * t75) * qJD(1) + m(3) * (-t93 * pkin(6) + t99) - t73 * mrSges(3,1) + t72 * mrSges(3,2) - t97;
t102 = t90 * t40 - t85 * t41;
t54 = -t65 * mrSges(4,1) + t66 * mrSges(4,2);
t26 = -t80 * pkin(3) - t79 * pkin(8) + t66 * t55 - t102;
t98 = t20 * mrSges(6,1) + t34 * t32 - m(6) * (-t30 * pkin(4) - t56 * pkin(9) + t58 * t52 + t26) - t21 * mrSges(6,2) - t35 * t33;
t94 = m(5) * t26 - t30 * mrSges(5,1) + t31 * mrSges(5,2) - t57 * t50 + t58 * t51 - t98;
t12 = m(4) * t102 + t80 * mrSges(4,1) - t48 * mrSges(4,3) - t66 * t54 + t81 * t59 - t94;
t7 = m(4) * t109 - t80 * mrSges(4,2) + t47 * mrSges(4,3) - t84 * t10 + t89 * t11 + t65 * t54 - t81 * t60;
t71 = (-mrSges(3,1) * t91 + mrSges(3,2) * t86) * qJD(1);
t4 = m(3) * (-t91 * g(3) - t111) - t72 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t71 * t108 + qJD(2) * t75 + t85 * t7 + t90 * t12;
t5 = m(3) * t104 - qJDD(2) * mrSges(3,2) + t73 * mrSges(3,3) - qJD(2) * t74 + t71 * t107 - t85 * t12 + t90 * t7;
t113 = t91 * t4 + t86 * t5;
t6 = m(2) * t105 + qJDD(1) * mrSges(2,1) - t93 * mrSges(2,2) - t114;
t1 = m(2) * t101 - t93 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t86 * t4 + t91 * t5;
t2 = [-m(1) * g(1) + t92 * t1 - t87 * t6, t1, t5, t7, t11, t14; -m(1) * g(2) + t87 * t1 + t92 * t6, t6, t4, t12, t10, t13; (-m(1) - m(2)) * g(3) + t113, -m(2) * g(3) + t113, t114, -t97, t94, -t98;];
f_new = t2;
