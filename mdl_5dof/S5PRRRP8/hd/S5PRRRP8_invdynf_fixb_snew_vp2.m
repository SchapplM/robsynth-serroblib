% Calculate vector of cutting forces with Newton-Euler
% S5PRRRP8
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 17:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRRRP8_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP8_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP8_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP8_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP8_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP8_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP8_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:58:20
% EndTime: 2019-12-05 16:58:22
% DurationCPUTime: 0.85s
% Computational Cost: add. (8445->132), mult. (15896->172), div. (0->0), fcn. (10297->10), ass. (0->72)
t63 = sin(pkin(9));
t65 = cos(pkin(9));
t53 = -t65 * g(1) - t63 * g(2);
t62 = -g(3) + qJDD(1);
t64 = sin(pkin(5));
t69 = sin(qJ(2));
t71 = cos(qJ(2));
t52 = t63 * g(1) - t65 * g(2);
t66 = cos(pkin(5));
t96 = t52 * t66;
t102 = -t69 * t53 + (t62 * t64 + t96) * t71;
t67 = sin(qJ(4));
t68 = sin(qJ(3));
t88 = qJD(2) * t68;
t97 = cos(qJ(4));
t47 = t67 * qJD(3) + t97 * t88;
t70 = cos(qJ(3));
t86 = qJD(2) * qJD(3);
t80 = t70 * t86;
t50 = t68 * qJDD(2) + t80;
t28 = t47 * qJD(4) - t97 * qJDD(3) + t67 * t50;
t87 = t70 * qJD(2);
t58 = qJD(4) - t87;
t36 = t58 * mrSges(5,1) - t47 * mrSges(5,3);
t81 = t68 * t86;
t51 = t70 * qJDD(2) - t81;
t43 = qJDD(4) - t51;
t46 = -t97 * qJD(3) + t67 * t88;
t31 = t46 * pkin(4) - t47 * qJ(5);
t37 = -t58 * mrSges(6,1) + t47 * mrSges(6,2);
t57 = t58 ^ 2;
t49 = (-pkin(3) * t70 - pkin(8) * t68) * qJD(2);
t72 = qJD(3) ^ 2;
t73 = qJD(2) ^ 2;
t95 = t64 * t69;
t83 = t71 * t53 + t62 * t95 + t69 * t96;
t25 = -t73 * pkin(2) + qJDD(2) * pkin(7) + t83;
t39 = -t64 * t52 + t66 * t62;
t91 = t70 * t25 + t68 * t39;
t19 = -t72 * pkin(3) + qJDD(3) * pkin(8) + t49 * t87 + t91;
t24 = -qJDD(2) * pkin(2) - t73 * pkin(7) - t102;
t21 = (-t50 - t80) * pkin(8) + (-t51 + t81) * pkin(3) + t24;
t92 = t97 * t19 + t67 * t21;
t85 = m(6) * (-t57 * pkin(4) + t43 * qJ(5) + 0.2e1 * qJD(5) * t58 - t46 * t31 + t92) + t58 * t37 + t43 * mrSges(6,3);
t32 = t46 * mrSges(6,1) - t47 * mrSges(6,3);
t90 = -t46 * mrSges(5,1) - t47 * mrSges(5,2) - t32;
t93 = -mrSges(5,3) - mrSges(6,2);
t11 = m(5) * t92 - t43 * mrSges(5,2) + t93 * t28 - t58 * t36 + t90 * t46 + t85;
t29 = -t46 * qJD(4) + t67 * qJDD(3) + t97 * t50;
t35 = -t58 * mrSges(5,2) - t46 * mrSges(5,3);
t38 = -t46 * mrSges(6,2) + t58 * mrSges(6,3);
t76 = -t67 * t19 + t97 * t21;
t99 = m(6) * (-t43 * pkin(4) - t57 * qJ(5) + t47 * t31 + qJDD(5) - t76);
t12 = m(5) * t76 - t99 + (t35 + t38) * t58 + t90 * t47 + (mrSges(5,1) + mrSges(6,1)) * t43 + t93 * t29;
t54 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t88;
t55 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t87;
t101 = m(4) * t24 - t51 * mrSges(4,1) + t50 * mrSges(4,2) + (t54 * t68 - t55 * t70) * qJD(2) + t67 * t11 + t97 * t12;
t79 = -t68 * t25 + t70 * t39;
t18 = -qJDD(3) * pkin(3) - t72 * pkin(8) + t49 * t88 - t79;
t84 = m(6) * (-0.2e1 * qJD(5) * t47 + (t46 * t58 - t29) * qJ(5) + (t47 * t58 + t28) * pkin(4) + t18) + t28 * mrSges(6,1) + t46 * t38;
t100 = m(5) * t18 + t28 * mrSges(5,1) + (mrSges(5,2) - mrSges(6,3)) * t29 + t46 * t35 + (t36 - t37) * t47 + t84;
t8 = m(3) * t102 + qJDD(2) * mrSges(3,1) - t73 * mrSges(3,2) - t101;
t98 = t71 * t8;
t48 = (-mrSges(4,1) * t70 + mrSges(4,2) * t68) * qJD(2);
t10 = m(4) * t79 + qJDD(3) * mrSges(4,1) - t50 * mrSges(4,3) + qJD(3) * t55 - t48 * t88 - t100;
t9 = m(4) * t91 - qJDD(3) * mrSges(4,2) + t51 * mrSges(4,3) - qJD(3) * t54 + t97 * t11 - t67 * t12 + t48 * t87;
t4 = m(3) * t83 - t73 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t68 * t10 + t70 * t9;
t6 = m(3) * t39 + t70 * t10 + t68 * t9;
t82 = m(2) * t62 + t4 * t95 + t66 * t6 + t64 * t98;
t2 = m(2) * t53 + t71 * t4 - t69 * t8;
t1 = m(2) * t52 - t64 * t6 + (t4 * t69 + t98) * t66;
t3 = [-m(1) * g(1) - t63 * t1 + t65 * t2, t2, t4, t9, t11, -t28 * mrSges(6,2) - t46 * t32 + t85; -m(1) * g(2) + t65 * t1 + t63 * t2, t1, t8, t10, t12, -t29 * mrSges(6,3) - t47 * t37 + t84; -m(1) * g(3) + t82, t82, t6, t101, t100, -t43 * mrSges(6,1) + t29 * mrSges(6,2) + t47 * t32 - t58 * t38 + t99;];
f_new = t3;
