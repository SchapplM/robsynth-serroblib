% Calculate vector of cutting forces with Newton-Euler
% S5RRRRR14
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
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
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
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRRR14_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR14_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR14_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR14_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR14_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR14_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR14_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 18:42:34
% EndTime: 2024-09-27 18:42:38
% DurationCPUTime: 2.05s
% Computational Cost: add. (53109->147), mult. (81218->209), div. (0->0), fcn. (57816->12), ass. (0->87)
t111 = -m(2) - m(3);
t75 = sin(pkin(5));
t110 = pkin(8) * t75;
t73 = qJD(1) + qJD(2);
t70 = t73 ^ 2;
t71 = qJDD(1) + qJDD(2);
t81 = sin(qJ(1));
t86 = cos(qJ(1));
t97 = t81 * g(1) - g(2) * t86;
t60 = qJDD(1) * pkin(1) + t97;
t87 = qJD(1) ^ 2;
t93 = -g(1) * t86 - g(2) * t81;
t62 = -pkin(1) * t87 + t93;
t80 = sin(qJ(2));
t85 = cos(qJ(2));
t94 = t85 * t60 - t62 * t80;
t41 = pkin(2) * t71 + t110 * t70 + t94;
t76 = cos(pkin(5));
t109 = t41 * t76;
t108 = t70 * t75 ^ 2;
t107 = t73 * t75;
t79 = sin(qJ(3));
t106 = t75 * t79;
t84 = cos(qJ(3));
t105 = t75 * t84;
t102 = qJD(3) * t73;
t54 = (t102 * t84 + t71 * t79) * t75;
t65 = t76 * t71 + qJDD(3);
t66 = t76 * t73 + qJD(3);
t103 = t80 * t60 + t85 * t62;
t42 = -pkin(2) * t70 + t110 * t71 + t103;
t95 = t84 * t109 - t79 * t42;
t24 = t65 * pkin(3) - t54 * pkin(9) + (pkin(3) * t79 * t108 + (pkin(9) * t66 * t73 - g(3)) * t75) * t84 + t95;
t100 = t84 ^ 2 * t108;
t99 = t73 * t106;
t53 = pkin(3) * t66 - pkin(9) * t99;
t55 = (-t102 * t79 + t71 * t84) * t75;
t91 = -g(3) * t106 + t79 * t109 + t84 * t42;
t25 = -pkin(3) * t100 + pkin(9) * t55 - t53 * t66 + t91;
t78 = sin(qJ(4));
t83 = cos(qJ(4));
t104 = t78 * t24 + t83 * t25;
t50 = mrSges(4,1) * t66 - mrSges(4,3) * t99;
t51 = mrSges(4,3) * t105 * t73 - mrSges(4,2) * t66;
t49 = (t78 * t84 + t79 * t83) * t107;
t30 = -qJD(4) * t49 - t54 * t78 + t55 * t83;
t48 = (-t78 * t79 + t83 * t84) * t107;
t31 = qJD(4) * t48 + t54 * t83 + t55 * t78;
t64 = qJD(4) + t66;
t43 = -mrSges(5,2) * t64 + mrSges(5,3) * t48;
t44 = mrSges(5,1) * t64 - mrSges(5,3) * t49;
t92 = -g(3) * t76 - t41 * t75;
t89 = -pkin(3) * t55 - pkin(9) * t100 + t53 * t99 + t92;
t77 = sin(qJ(5));
t82 = cos(qJ(5));
t36 = t48 * t77 + t49 * t82;
t19 = -qJD(5) * t36 + t30 * t82 - t31 * t77;
t35 = t48 * t82 - t49 * t77;
t20 = qJD(5) * t35 + t30 * t77 + t31 * t82;
t61 = qJD(5) + t64;
t32 = -mrSges(6,2) * t61 + mrSges(6,3) * t35;
t33 = mrSges(6,1) * t61 - mrSges(6,3) * t36;
t45 = pkin(4) * t64 - pkin(10) * t49;
t47 = t48 ^ 2;
t90 = -t19 * mrSges(6,1) - t35 * t32 + m(6) * (-pkin(4) * t30 - pkin(10) * t47 + t45 * t49 + t89) + t20 * mrSges(6,2) + t36 * t33;
t88 = m(5) * t89 - t30 * mrSges(5,1) + t31 * mrSges(5,2) - t48 * t43 + t49 * t44 + t90;
t12 = (t50 * t79 - t51 * t84) * t107 + m(4) * t92 + t54 * mrSges(4,2) - t55 * mrSges(4,1) + t88;
t63 = qJDD(4) + t65;
t96 = t83 * t24 - t78 * t25;
t15 = (t48 * t64 - t31) * pkin(10) + (t48 * t49 + t63) * pkin(4) + t96;
t16 = -pkin(4) * t47 + pkin(10) * t30 - t45 * t64 + t104;
t28 = -mrSges(6,1) * t35 + mrSges(6,2) * t36;
t59 = qJDD(5) + t63;
t13 = m(6) * (t15 * t82 - t16 * t77) - t20 * mrSges(6,3) + t59 * mrSges(6,1) - t36 * t28 + t61 * t32;
t14 = m(6) * (t15 * t77 + t16 * t82) + t19 * mrSges(6,3) - t59 * mrSges(6,2) + t35 * t28 - t61 * t33;
t37 = -mrSges(5,1) * t48 + mrSges(5,2) * t49;
t10 = m(5) * t104 - t63 * mrSges(5,2) + t30 * mrSges(5,3) - t77 * t13 + t82 * t14 + t48 * t37 - t64 * t44;
t9 = m(5) * t96 + t63 * mrSges(5,1) - t31 * mrSges(5,3) + t82 * t13 + t77 * t14 - t49 * t37 + t64 * t43;
t98 = (-mrSges(4,1) * t84 + mrSges(4,2) * t79) * t107 ^ 2;
t7 = m(4) * (-g(3) * t105 + t95) - t54 * mrSges(4,3) + t65 * mrSges(4,1) - t79 * t98 + t66 * t51 + t78 * t10 + t83 * t9;
t8 = m(4) * t91 - t65 * mrSges(4,2) + t55 * mrSges(4,3) + t83 * t10 - t66 * t50 - t78 * t9 + t84 * t98;
t101 = t7 * t105 + t8 * t106 + t76 * t12;
t4 = m(3) * t103 - t70 * mrSges(3,1) - t71 * mrSges(3,2) - t79 * t7 + t84 * t8;
t3 = m(3) * t94 + t71 * mrSges(3,1) - t70 * mrSges(3,2) - t75 * t12 + (t7 * t84 + t79 * t8) * t76;
t2 = m(2) * t93 - t87 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t80 * t3 + t85 * t4;
t1 = m(2) * t97 + qJDD(1) * mrSges(2,1) - t87 * mrSges(2,2) + t85 * t3 + t80 * t4;
t5 = [-m(1) * g(1) - t1 * t81 + t2 * t86, t2, t4, t8, t10, t14; -m(1) * g(2) + t1 * t86 + t2 * t81, t1, t3, t7, t9, t13; (-m(1) + t111) * g(3) + t101, g(3) * t111 + t101, -m(3) * g(3) + t101, t12, t88, t90;];
f_new = t5;
