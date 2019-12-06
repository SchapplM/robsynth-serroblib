% Calculate vector of cutting forces with Newton-Euler
% S5RPPRR2
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPPRR2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR2_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR2_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR2_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR2_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR2_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:39:37
% EndTime: 2019-12-05 17:39:39
% DurationCPUTime: 0.90s
% Computational Cost: add. (8383->126), mult. (18803->159), div. (0->0), fcn. (12583->8), ass. (0->71)
t70 = qJD(1) ^ 2;
t66 = sin(qJ(1));
t69 = cos(qJ(1));
t87 = t66 * g(1) - t69 * g(2);
t77 = -t70 * qJ(2) + qJDD(2) - t87;
t96 = -pkin(1) - qJ(3);
t103 = -(2 * qJD(1) * qJD(3)) + t96 * qJDD(1) + t77;
t62 = sin(pkin(8));
t58 = t62 ^ 2;
t63 = cos(pkin(8));
t94 = t63 ^ 2 + t58;
t86 = t94 * mrSges(4,3);
t83 = -t69 * g(1) - t66 * g(2);
t102 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t83;
t101 = -m(2) - m(3);
t100 = pkin(3) * t70;
t99 = mrSges(4,2) * t63;
t98 = mrSges(2,1) - mrSges(3,2);
t97 = -mrSges(2,2) + mrSges(3,3);
t90 = t62 * g(3) + t103 * t63;
t22 = (-pkin(6) * qJDD(1) - t62 * t100) * t63 + t90;
t84 = -t63 * g(3) + t103 * t62;
t92 = qJDD(1) * t62;
t23 = -pkin(6) * t92 - t58 * t100 + t84;
t65 = sin(qJ(4));
t68 = cos(qJ(4));
t95 = t65 * t22 + t68 * t23;
t82 = -t62 * t68 - t63 * t65;
t46 = t82 * qJD(1);
t93 = t46 * qJD(4);
t81 = -t62 * t65 + t63 * t68;
t47 = t81 * qJD(1);
t31 = -t46 * mrSges(5,1) + t47 * mrSges(5,2);
t34 = t81 * qJDD(1) + t93;
t39 = -qJD(4) * mrSges(5,2) + t46 * mrSges(5,3);
t64 = sin(qJ(5));
t67 = cos(qJ(5));
t85 = t68 * t22 - t65 * t23;
t10 = (-t34 + t93) * pkin(7) + (t46 * t47 + qJDD(4)) * pkin(4) + t85;
t33 = -t47 * qJD(4) + t82 * qJDD(1);
t41 = qJD(4) * pkin(4) - t47 * pkin(7);
t45 = t46 ^ 2;
t11 = -t45 * pkin(4) + t33 * pkin(7) - qJD(4) * t41 + t95;
t28 = t67 * t46 - t64 * t47;
t16 = t28 * qJD(5) + t64 * t33 + t67 * t34;
t29 = t64 * t46 + t67 * t47;
t18 = -t28 * mrSges(6,1) + t29 * mrSges(6,2);
t60 = qJD(4) + qJD(5);
t24 = -t60 * mrSges(6,2) + t28 * mrSges(6,3);
t57 = qJDD(4) + qJDD(5);
t8 = m(6) * (t67 * t10 - t64 * t11) - t16 * mrSges(6,3) + t57 * mrSges(6,1) - t29 * t18 + t60 * t24;
t15 = -t29 * qJD(5) + t67 * t33 - t64 * t34;
t25 = t60 * mrSges(6,1) - t29 * mrSges(6,3);
t9 = m(6) * (t64 * t10 + t67 * t11) + t15 * mrSges(6,3) - t57 * mrSges(6,2) + t28 * t18 - t60 * t25;
t5 = m(5) * t85 + qJDD(4) * mrSges(5,1) - t34 * mrSges(5,3) + qJD(4) * t39 - t47 * t31 + t64 * t9 + t67 * t8;
t40 = qJD(4) * mrSges(5,1) - t47 * mrSges(5,3);
t6 = m(5) * t95 - qJDD(4) * mrSges(5,2) + t33 * mrSges(5,3) - qJD(4) * t40 + t46 * t31 - t64 * t8 + t67 * t9;
t79 = -qJDD(1) * mrSges(4,3) - t70 * (mrSges(4,1) * t62 + t99);
t3 = m(4) * t90 + t68 * t5 + t65 * t6 + t79 * t63;
t4 = m(4) * t84 - t65 * t5 + t68 * t6 + t79 * t62;
t88 = -t62 * t3 + t63 * t4;
t78 = -m(3) * (-qJDD(1) * pkin(1) + t77) - t63 * t3 - t62 * t4;
t75 = qJDD(3) + t102;
t73 = pkin(3) * t92 + (-t94 * pkin(6) + t96) * t70 + t75;
t76 = t15 * mrSges(6,1) + t28 * t24 - m(6) * (-t33 * pkin(4) - t45 * pkin(7) + t47 * t41 + t73) - t16 * mrSges(6,2) - t29 * t25;
t74 = -m(5) * t73 + t33 * mrSges(5,1) - t34 * mrSges(5,2) + t46 * t39 - t47 * t40 + t76;
t72 = m(4) * (t96 * t70 + t75) + mrSges(4,1) * t92 + qJDD(1) * t99 - t74;
t71 = m(3) * (t70 * pkin(1) - t102) - t72;
t7 = -t71 + (-t86 - t98) * t70 + t97 * qJDD(1) + m(2) * t83;
t1 = m(2) * t87 + t98 * qJDD(1) + t97 * t70 + t78;
t2 = [-m(1) * g(1) - t66 * t1 + t69 * t7, t7, -m(3) * g(3) + t88, t4, t6, t9; -m(1) * g(2) + t69 * t1 + t66 * t7, t1, -qJDD(1) * mrSges(3,3) + t71 + (-mrSges(3,2) + t86) * t70, t3, t5, t8; (-m(1) + t101) * g(3) + t88, t101 * g(3) + t88, qJDD(1) * mrSges(3,2) - t70 * mrSges(3,3) - t78, -t70 * t86 + t72, -t74, -t76;];
f_new = t2;
