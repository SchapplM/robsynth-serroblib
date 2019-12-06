% Calculate vector of cutting forces with Newton-Euler
% S5RPRPR4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-05 17:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRPR4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR4_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR4_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR4_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR4_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:53:16
% EndTime: 2019-12-05 17:53:19
% DurationCPUTime: 1.33s
% Computational Cost: add. (14441->138), mult. (31070->188), div. (0->0), fcn. (19687->10), ass. (0->73)
t76 = sin(qJ(3));
t79 = cos(qJ(3));
t95 = qJD(1) * qJD(3);
t91 = t79 * t95;
t56 = t76 * qJDD(1) + t91;
t57 = t79 * qJDD(1) - t76 * t95;
t97 = qJD(1) * t76;
t59 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t97;
t96 = qJD(1) * t79;
t60 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t96;
t81 = qJD(1) ^ 2;
t71 = sin(pkin(9));
t73 = cos(pkin(9));
t39 = -t71 * t56 + t73 * t57;
t40 = t73 * t56 + t71 * t57;
t47 = (-t71 * t76 + t73 * t79) * qJD(1);
t41 = -qJD(3) * mrSges(5,2) + t47 * mrSges(5,3);
t48 = (t71 * t79 + t73 * t76) * qJD(1);
t42 = qJD(3) * mrSges(5,1) - t48 * mrSges(5,3);
t58 = qJD(3) * pkin(3) - qJ(4) * t97;
t69 = t79 ^ 2;
t77 = sin(qJ(1));
t80 = cos(qJ(1));
t98 = t80 * g(2) + t77 * g(3);
t53 = qJDD(1) * pkin(1) + t98;
t92 = t77 * g(2) - t80 * g(3);
t55 = -t81 * pkin(1) + t92;
t72 = sin(pkin(8));
t74 = cos(pkin(8));
t89 = t74 * t53 - t72 * t55;
t86 = -qJDD(1) * pkin(2) - t89;
t83 = -t57 * pkin(3) + qJDD(4) + t58 * t97 + (-qJ(4) * t69 - pkin(6)) * t81 + t86;
t75 = sin(qJ(5));
t78 = cos(qJ(5));
t33 = t75 * t47 + t78 * t48;
t18 = -t33 * qJD(5) + t78 * t39 - t75 * t40;
t32 = t78 * t47 - t75 * t48;
t19 = t32 * qJD(5) + t75 * t39 + t78 * t40;
t68 = qJD(3) + qJD(5);
t29 = -t68 * mrSges(6,2) + t32 * mrSges(6,3);
t30 = t68 * mrSges(6,1) - t33 * mrSges(6,3);
t43 = qJD(3) * pkin(4) - t48 * pkin(7);
t46 = t47 ^ 2;
t85 = t18 * mrSges(6,1) + t32 * t29 - m(6) * (-t39 * pkin(4) - t46 * pkin(7) + t48 * t43 + t83) - t19 * mrSges(6,2) - t33 * t30;
t84 = -m(5) * t83 + t39 * mrSges(5,1) - t40 * mrSges(5,2) + t47 * t41 - t48 * t42 + t85;
t101 = (t76 * t59 - t79 * t60) * qJD(1) + m(4) * (-t81 * pkin(6) + t86) - t57 * mrSges(4,1) + t56 * mrSges(4,2) - t84;
t99 = t72 * t53 + t74 * t55;
t35 = -t81 * pkin(2) + qJDD(1) * pkin(6) + t99;
t70 = -g(1) + qJDD(2);
t100 = t79 * t35 + t76 * t70;
t54 = (-mrSges(4,1) * t79 + mrSges(4,2) * t76) * qJD(1);
t90 = -t76 * t35 + t79 * t70;
t25 = (-t56 + t91) * qJ(4) + (t76 * t79 * t81 + qJDD(3)) * pkin(3) + t90;
t26 = -t69 * t81 * pkin(3) + t57 * qJ(4) - qJD(3) * t58 + t100;
t88 = -0.2e1 * qJD(4) * t48 + t73 * t25 - t71 * t26;
t13 = (qJD(3) * t47 - t40) * pkin(7) + (t47 * t48 + qJDD(3)) * pkin(4) + t88;
t93 = 0.2e1 * qJD(4) * t47 + t71 * t25 + t73 * t26;
t14 = -t46 * pkin(4) + t39 * pkin(7) - qJD(3) * t43 + t93;
t24 = -mrSges(6,1) * t32 + mrSges(6,2) * t33;
t67 = qJDD(3) + qJDD(5);
t11 = m(6) * (t78 * t13 - t75 * t14) - t19 * mrSges(6,3) + t67 * mrSges(6,1) - t33 * t24 + t68 * t29;
t12 = m(6) * (t75 * t13 + t78 * t14) + t18 * mrSges(6,3) - t67 * mrSges(6,2) + t32 * t24 - t68 * t30;
t37 = -t47 * mrSges(5,1) + t48 * mrSges(5,2);
t8 = m(5) * t88 + qJDD(3) * mrSges(5,1) - t40 * mrSges(5,3) + qJD(3) * t41 + t78 * t11 + t75 * t12 - t48 * t37;
t9 = m(5) * t93 - qJDD(3) * mrSges(5,2) + t39 * mrSges(5,3) - qJD(3) * t42 - t75 * t11 + t78 * t12 + t47 * t37;
t6 = m(4) * t90 + qJDD(3) * mrSges(4,1) - t56 * mrSges(4,3) + qJD(3) * t60 - t54 * t97 + t71 * t9 + t73 * t8;
t7 = m(4) * t100 - qJDD(3) * mrSges(4,2) + t57 * mrSges(4,3) - qJD(3) * t59 + t54 * t96 - t71 * t8 + t73 * t9;
t94 = m(3) * t70 + t79 * t6 + t76 * t7;
t10 = m(3) * t89 + qJDD(1) * mrSges(3,1) - t81 * mrSges(3,2) - t101;
t3 = m(3) * t99 - t81 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t76 * t6 + t79 * t7;
t2 = m(2) * t98 + qJDD(1) * mrSges(2,1) - t81 * mrSges(2,2) + t74 * t10 + t72 * t3;
t1 = m(2) * t92 - t81 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t72 * t10 + t74 * t3;
t4 = [(-m(1) - m(2)) * g(1) + t94, t1, t3, t7, t9, t12; -m(1) * g(2) - t77 * t1 - t80 * t2, t2, t10, t6, t8, t11; -m(1) * g(3) + t80 * t1 - t77 * t2, -m(2) * g(1) + t94, t94, t101, -t84, -t85;];
f_new = t4;
