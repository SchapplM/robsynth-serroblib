% Calculate vector of cutting forces with Newton-Euler
% S5PRRRR11
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
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRRRR11_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR11_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR11_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR11_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR11_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR11_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR11_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 21:45:22
% EndTime: 2024-09-27 21:45:24
% DurationCPUTime: 1.00s
% Computational Cost: add. (30887->137), mult. (76346->199), div. (0->0), fcn. (57816->12), ass. (0->84)
t78 = sin(pkin(5));
t112 = pkin(7) * t78;
t89 = qJD(2) ^ 2;
t77 = sin(pkin(10));
t79 = cos(pkin(10));
t61 = t77 * g(1) - t79 * g(2);
t62 = -t79 * g(1) - t77 * g(2);
t84 = sin(qJ(2));
t88 = cos(qJ(2));
t95 = t88 * t61 - t84 * t62;
t41 = qJDD(2) * pkin(2) + t89 * t112 + t95;
t80 = cos(pkin(5));
t111 = t41 * t80;
t110 = t78 ^ 2 * t89;
t83 = sin(qJ(3));
t109 = t78 * t83;
t87 = cos(qJ(3));
t108 = t78 * t87;
t105 = qJD(2) * t78;
t100 = t87 * t105;
t104 = qJD(2) * qJD(3);
t54 = (qJDD(2) * t83 + t87 * t104) * t78;
t67 = t80 * qJDD(2) + qJDD(3);
t68 = t80 * qJD(2) + qJD(3);
t106 = t84 * t61 + t88 * t62;
t42 = -t89 * pkin(2) + qJDD(2) * t112 + t106;
t76 = -g(3) + qJDD(1);
t93 = t76 * t108 + t87 * t111 - t83 * t42;
t24 = (t68 * t100 - t54) * pkin(8) + (t83 * t87 * t110 + t67) * pkin(3) + t93;
t102 = t76 * t109 + t83 * t111 + t87 * t42;
t103 = t87 ^ 2 * t110;
t101 = t83 * t105;
t53 = t68 * pkin(3) - pkin(8) * t101;
t55 = (qJDD(2) * t87 - t83 * t104) * t78;
t25 = -pkin(3) * t103 + t55 * pkin(8) - t68 * t53 + t102;
t82 = sin(qJ(4));
t86 = cos(qJ(4));
t107 = t82 * t24 + t86 * t25;
t99 = (-mrSges(4,1) * t87 + mrSges(4,2) * t83) * t105 ^ 2;
t66 = qJD(4) + t68;
t50 = t68 * mrSges(4,1) - mrSges(4,3) * t101;
t51 = -t68 * mrSges(4,2) + mrSges(4,3) * t100;
t49 = (t82 * t87 + t83 * t86) * t105;
t30 = -t49 * qJD(4) - t82 * t54 + t86 * t55;
t48 = (-t82 * t83 + t86 * t87) * t105;
t31 = t48 * qJD(4) + t86 * t54 + t82 * t55;
t43 = -t66 * mrSges(5,2) + t48 * mrSges(5,3);
t44 = t66 * mrSges(5,1) - t49 * mrSges(5,3);
t96 = -t78 * t41 + t80 * t76;
t91 = -t55 * pkin(3) - pkin(8) * t103 + t53 * t101 + t96;
t81 = sin(qJ(5));
t85 = cos(qJ(5));
t36 = t81 * t48 + t85 * t49;
t19 = -t36 * qJD(5) + t85 * t30 - t81 * t31;
t35 = t85 * t48 - t81 * t49;
t20 = t35 * qJD(5) + t81 * t30 + t85 * t31;
t60 = qJD(5) + t66;
t32 = -t60 * mrSges(6,2) + t35 * mrSges(6,3);
t33 = t60 * mrSges(6,1) - t36 * mrSges(6,3);
t45 = t66 * pkin(4) - t49 * pkin(9);
t47 = t48 ^ 2;
t92 = -t19 * mrSges(6,1) - t35 * t32 + m(6) * (-t30 * pkin(4) - t47 * pkin(9) + t49 * t45 + t91) + t20 * mrSges(6,2) + t36 * t33;
t90 = m(5) * t91 - t30 * mrSges(5,1) + t31 * mrSges(5,2) - t48 * t43 + t49 * t44 + t92;
t12 = (t50 * t83 - t51 * t87) * t105 + m(4) * t96 + t54 * mrSges(4,2) - t55 * mrSges(4,1) + t90;
t65 = qJDD(4) + t67;
t97 = t86 * t24 - t82 * t25;
t15 = (t48 * t66 - t31) * pkin(9) + (t48 * t49 + t65) * pkin(4) + t97;
t16 = -t47 * pkin(4) + t30 * pkin(9) - t66 * t45 + t107;
t27 = -t35 * mrSges(6,1) + t36 * mrSges(6,2);
t59 = qJDD(5) + t65;
t13 = m(6) * (t85 * t15 - t81 * t16) - t20 * mrSges(6,3) + t59 * mrSges(6,1) - t36 * t27 + t60 * t32;
t14 = m(6) * (t81 * t15 + t85 * t16) + t19 * mrSges(6,3) - t59 * mrSges(6,2) + t35 * t27 - t60 * t33;
t37 = -t48 * mrSges(5,1) + t49 * mrSges(5,2);
t10 = m(5) * t107 - t65 * mrSges(5,2) + t30 * mrSges(5,3) - t81 * t13 + t85 * t14 + t48 * t37 - t66 * t44;
t9 = m(5) * t97 + t65 * mrSges(5,1) - t31 * mrSges(5,3) + t85 * t13 + t81 * t14 - t49 * t37 + t66 * t43;
t7 = m(4) * t93 + t67 * mrSges(4,1) - t54 * mrSges(4,3) + t82 * t10 + t68 * t51 - t83 * t99 + t86 * t9;
t8 = m(4) * t102 - t67 * mrSges(4,2) + t55 * mrSges(4,3) + t86 * t10 - t68 * t50 - t82 * t9 + t87 * t99;
t98 = m(3) * t76 + t7 * t108 + t8 * t109 + t80 * t12;
t94 = m(2) * t76 + t98;
t4 = m(3) * t106 - t89 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t83 * t7 + t87 * t8;
t3 = m(3) * t95 + qJDD(2) * mrSges(3,1) - t89 * mrSges(3,2) - t78 * t12 + (t87 * t7 + t83 * t8) * t80;
t2 = m(2) * t62 - t84 * t3 + t88 * t4;
t1 = m(2) * t61 + t88 * t3 + t84 * t4;
t5 = [-m(1) * g(1) - t77 * t1 + t79 * t2, t2, t4, t8, t10, t14; -m(1) * g(2) + t79 * t1 + t77 * t2, t1, t3, t7, t9, t13; -m(1) * g(3) + t94, t94, t98, t12, t90, t92;];
f_new = t5;
