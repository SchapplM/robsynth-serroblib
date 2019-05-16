% Calculate vector of cutting forces with Newton-Euler
% S6PPRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-05-04 20:37
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PPRRRP2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP2_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRP2_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRP2_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRP2_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRP2_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:32:03
% EndTime: 2019-05-04 20:32:07
% DurationCPUTime: 2.34s
% Computational Cost: add. (33578->142), mult. (59974->188), div. (0->0), fcn. (45809->14), ass. (0->83)
t70 = sin(pkin(11));
t74 = cos(pkin(11));
t60 = -t74 * g(1) - t70 * g(2);
t69 = sin(pkin(12));
t73 = cos(pkin(12));
t59 = t70 * g(1) - t74 * g(2);
t68 = -g(3) + qJDD(1);
t72 = sin(pkin(6));
t76 = cos(pkin(6));
t88 = t59 * t76 + t68 * t72;
t34 = -t69 * t60 + t88 * t73;
t49 = -t72 * t59 + t76 * t68 + qJDD(2);
t71 = sin(pkin(7));
t75 = cos(pkin(7));
t117 = t34 * t75 + t49 * t71;
t83 = qJD(3) ^ 2;
t35 = t73 * t60 + t88 * t69;
t79 = sin(qJ(3));
t81 = cos(qJ(3));
t95 = t117 * t79 + t81 * t35;
t28 = -t83 * pkin(3) + qJDD(3) * pkin(9) + t95;
t30 = -t71 * t34 + t75 * t49;
t78 = sin(qJ(4));
t80 = cos(qJ(4));
t103 = t80 * t28 + t78 * t30;
t109 = cos(qJ(5));
t100 = qJD(3) * t78;
t77 = sin(qJ(5));
t53 = -t109 * qJD(4) + t77 * t100;
t54 = t77 * qJD(4) + t109 * t100;
t42 = t53 * mrSges(7,1) - t54 * mrSges(7,3);
t102 = -t53 * mrSges(6,1) - t54 * mrSges(6,2) - t42;
t56 = (-pkin(4) * t80 - pkin(10) * t78) * qJD(3);
t82 = qJD(4) ^ 2;
t99 = t80 * qJD(3);
t22 = -t82 * pkin(4) + qJDD(4) * pkin(10) + t56 * t99 + t103;
t115 = t117 * t81 - t79 * t35;
t27 = -qJDD(3) * pkin(3) - t83 * pkin(9) - t115;
t98 = qJD(3) * qJD(4);
t92 = t80 * t98;
t57 = t78 * qJDD(3) + t92;
t93 = t78 * t98;
t58 = t80 * qJDD(3) - t93;
t24 = (-t57 - t92) * pkin(10) + (-t58 + t93) * pkin(4) + t27;
t104 = t109 * t22 + t77 * t24;
t105 = -mrSges(6,3) - mrSges(7,2);
t38 = t54 * qJD(5) - t109 * qJDD(4) + t77 * t57;
t64 = qJD(5) - t99;
t46 = t64 * mrSges(6,1) - t54 * mrSges(6,3);
t52 = qJDD(5) - t58;
t41 = t53 * pkin(5) - t54 * qJ(6);
t47 = -t64 * mrSges(7,1) + t54 * mrSges(7,2);
t63 = t64 ^ 2;
t97 = m(7) * (-t63 * pkin(5) + t52 * qJ(6) + 0.2e1 * qJD(6) * t64 - t53 * t41 + t104) + t64 * t47 + t52 * mrSges(7,3);
t14 = m(6) * t104 - t52 * mrSges(6,2) + t102 * t53 + t105 * t38 - t64 * t46 + t97;
t86 = t109 * t24 - t77 * t22;
t112 = m(7) * (-t52 * pkin(5) - t63 * qJ(6) + t54 * t41 + qJDD(6) - t86);
t39 = -t53 * qJD(5) + t77 * qJDD(4) + t109 * t57;
t45 = -t64 * mrSges(6,2) - t53 * mrSges(6,3);
t48 = -t53 * mrSges(7,2) + t64 * mrSges(7,3);
t15 = m(6) * t86 - t112 + (t45 + t48) * t64 + t102 * t54 + (mrSges(6,1) + mrSges(7,1)) * t52 + t105 * t39;
t55 = (-mrSges(5,1) * t80 + mrSges(5,2) * t78) * qJD(3);
t61 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t100;
t12 = m(5) * t103 - qJDD(4) * mrSges(5,2) + t58 * mrSges(5,3) - qJD(4) * t61 + t109 * t14 - t77 * t15 + t55 * t99;
t91 = -t78 * t28 + t80 * t30;
t21 = -qJDD(4) * pkin(4) - t82 * pkin(10) + t56 * t100 - t91;
t96 = m(7) * (-0.2e1 * qJD(6) * t54 + (t53 * t64 - t39) * qJ(6) + (t54 * t64 + t38) * pkin(5) + t21) + t38 * mrSges(7,1) + t53 * t48;
t113 = m(6) * t21 + t38 * mrSges(6,1) + (t46 - t47) * t54 + (mrSges(6,2) - mrSges(7,3)) * t39 + t53 * t45 + t96;
t62 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t99;
t13 = m(5) * t91 + qJDD(4) * mrSges(5,1) - t57 * mrSges(5,3) + qJD(4) * t62 - t55 * t100 - t113;
t10 = m(4) * t30 + t78 * t12 + t80 * t13;
t114 = m(5) * t27 - t58 * mrSges(5,1) + t57 * mrSges(5,2) + (t61 * t78 - t62 * t80) * qJD(3) + t109 * t15 + t77 * t14;
t11 = m(4) * t115 + qJDD(3) * mrSges(4,1) - t83 * mrSges(4,2) - t114;
t9 = m(4) * t95 - t83 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t80 * t12 - t78 * t13;
t90 = t11 * t81 + t79 * t9;
t4 = m(3) * t34 - t71 * t10 + t90 * t75;
t8 = m(3) * t35 - t79 * t11 + t81 * t9;
t116 = t4 * t73 + t69 * t8;
t6 = m(3) * t49 + t75 * t10 + t90 * t71;
t94 = m(2) * t68 + t116 * t72 + t76 * t6;
t2 = m(2) * t60 - t69 * t4 + t73 * t8;
t1 = m(2) * t59 + t116 * t76 - t72 * t6;
t3 = [-m(1) * g(1) - t70 * t1 + t74 * t2, t2, t8, t9, t12, t14, -t38 * mrSges(7,2) - t53 * t42 + t97; -m(1) * g(2) + t74 * t1 + t70 * t2, t1, t4, t11, t13, t15, -t39 * mrSges(7,3) - t54 * t47 + t96; -m(1) * g(3) + t94, t94, t6, t10, t114, t113, -t52 * mrSges(7,1) + t39 * mrSges(7,2) + t54 * t42 - t64 * t48 + t112;];
f_new  = t3;
