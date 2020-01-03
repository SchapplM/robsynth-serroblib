% Calculate vector of cutting forces with Newton-Euler
% S5RPRRR10
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRRR10_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR10_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR10_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR10_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR10_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR10_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR10_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR10_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR10_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:09:48
% EndTime: 2019-12-31 19:09:52
% DurationCPUTime: 1.84s
% Computational Cost: add. (20783->155), mult. (49096->202), div. (0->0), fcn. (35932->10), ass. (0->86)
t88 = qJD(1) ^ 2;
t78 = cos(pkin(9));
t76 = t78 ^ 2;
t77 = sin(pkin(9));
t107 = t77 ^ 2 + t76;
t113 = t107 * mrSges(3,3);
t104 = qJD(1) * qJD(2);
t102 = -g(3) * t78 - 0.2e1 * t77 * t104;
t111 = pkin(2) * t88;
t82 = sin(qJ(1));
t86 = cos(qJ(1));
t98 = -g(1) * t86 - g(2) * t82;
t68 = -pkin(1) * t88 + qJDD(1) * qJ(2) + t98;
t44 = (-pkin(6) * qJDD(1) + t78 * t111 - t68) * t77 + t102;
t105 = qJDD(1) * t78;
t99 = -g(3) * t77 + (0.2e1 * t104 + t68) * t78;
t45 = pkin(6) * t105 - t76 * t111 + t99;
t81 = sin(qJ(3));
t85 = cos(qJ(3));
t100 = t44 * t85 - t81 * t45;
t110 = t77 * t81;
t66 = (t78 * t85 - t110) * qJD(1);
t95 = t77 * t85 + t78 * t81;
t67 = t95 * qJD(1);
t49 = -mrSges(4,1) * t66 + mrSges(4,2) * t67;
t106 = qJD(3) * t66;
t55 = t95 * qJDD(1) + t106;
t59 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t66;
t52 = -pkin(3) * t66 - pkin(7) * t67;
t87 = qJD(3) ^ 2;
t23 = -qJDD(3) * pkin(3) - pkin(7) * t87 + t67 * t52 - t100;
t80 = sin(qJ(4));
t84 = cos(qJ(4));
t58 = qJD(3) * t80 + t67 * t84;
t32 = -qJD(4) * t58 + qJDD(3) * t84 - t55 * t80;
t57 = qJD(3) * t84 - t67 * t80;
t33 = qJD(4) * t57 + qJDD(3) * t80 + t55 * t84;
t64 = qJD(4) - t66;
t39 = -mrSges(5,2) * t64 + mrSges(5,3) * t57;
t40 = mrSges(5,1) * t64 - mrSges(5,3) * t58;
t79 = sin(qJ(5));
t83 = cos(qJ(5));
t35 = t57 * t79 + t58 * t83;
t20 = -qJD(5) * t35 + t32 * t83 - t33 * t79;
t34 = t57 * t83 - t58 * t79;
t21 = qJD(5) * t34 + t32 * t79 + t33 * t83;
t62 = qJD(5) + t64;
t30 = -mrSges(6,2) * t62 + mrSges(6,3) * t34;
t31 = mrSges(6,1) * t62 - mrSges(6,3) * t35;
t43 = pkin(4) * t64 - pkin(8) * t58;
t56 = t57 ^ 2;
t93 = t20 * mrSges(6,1) + t34 * t30 - m(6) * (-pkin(4) * t32 - pkin(8) * t56 + t43 * t58 + t23) - t21 * mrSges(6,2) - t35 * t31;
t89 = m(5) * t23 - t32 * mrSges(5,1) + t33 * mrSges(5,2) - t57 * t39 + t58 * t40 - t93;
t12 = m(4) * t100 + qJDD(3) * mrSges(4,1) - t55 * mrSges(4,3) + qJD(3) * t59 - t67 * t49 - t89;
t108 = t81 * t44 + t85 * t45;
t24 = -pkin(3) * t87 + qJDD(3) * pkin(7) + t52 * t66 + t108;
t63 = t67 * qJD(3);
t54 = -qJDD(1) * t110 + t85 * t105 - t63;
t103 = g(1) * t82 - t86 * g(2);
t97 = qJDD(2) - t103;
t90 = (-pkin(2) * t78 - pkin(1)) * qJDD(1) + (-t107 * pkin(6) - qJ(2)) * t88 + t97;
t27 = (-t55 - t106) * pkin(7) + (-t54 + t63) * pkin(3) + t90;
t101 = -t24 * t80 + t84 * t27;
t51 = qJDD(4) - t54;
t15 = (t57 * t64 - t33) * pkin(8) + (t57 * t58 + t51) * pkin(4) + t101;
t109 = t84 * t24 + t80 * t27;
t16 = -pkin(4) * t56 + pkin(8) * t32 - t43 * t64 + t109;
t29 = -mrSges(6,1) * t34 + mrSges(6,2) * t35;
t47 = qJDD(5) + t51;
t13 = m(6) * (t15 * t83 - t16 * t79) - t21 * mrSges(6,3) + t47 * mrSges(6,1) - t35 * t29 + t62 * t30;
t14 = m(6) * (t15 * t79 + t16 * t83) + t20 * mrSges(6,3) - t47 * mrSges(6,2) + t34 * t29 - t62 * t31;
t36 = -mrSges(5,1) * t57 + mrSges(5,2) * t58;
t10 = m(5) * t101 + t51 * mrSges(5,1) - t33 * mrSges(5,3) + t83 * t13 + t79 * t14 - t58 * t36 + t64 * t39;
t11 = m(5) * t109 - t51 * mrSges(5,2) + t32 * mrSges(5,3) - t79 * t13 + t83 * t14 + t57 * t36 - t64 * t40;
t60 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t67;
t7 = m(4) * t108 - qJDD(3) * mrSges(4,2) + t54 * mrSges(4,3) - qJD(3) * t60 - t80 * t10 + t84 * t11 + t66 * t49;
t96 = -mrSges(3,1) * t78 + mrSges(3,2) * t77;
t94 = mrSges(3,3) * qJDD(1) + t88 * t96;
t4 = m(3) * t102 + t81 * t7 + t85 * t12 + (-m(3) * t68 - t94) * t77;
t5 = m(3) * t99 - t81 * t12 + t85 * t7 + t94 * t78;
t112 = t78 * t4 + t77 * t5;
t92 = -m(4) * t90 + t54 * mrSges(4,1) - t55 * mrSges(4,2) - t84 * t10 - t80 * t11 + t66 * t59 - t67 * t60;
t91 = m(3) * (-qJDD(1) * pkin(1) - qJ(2) * t88 + t97) - t92;
t6 = m(2) * t103 + (-mrSges(2,2) + t113) * t88 + (mrSges(2,1) - t96) * qJDD(1) - t91;
t1 = m(2) * t98 - t88 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t77 * t4 + t78 * t5;
t2 = [-m(1) * g(1) + t1 * t86 - t6 * t82, t1, t5, t7, t11, t14; -m(1) * g(2) + t1 * t82 + t6 * t86, t6, t4, t12, t10, t13; (-m(1) - m(2)) * g(3) + t112, -m(2) * g(3) + t112, t96 * qJDD(1) - t88 * t113 + t91, -t92, t89, -t93;];
f_new = t2;
