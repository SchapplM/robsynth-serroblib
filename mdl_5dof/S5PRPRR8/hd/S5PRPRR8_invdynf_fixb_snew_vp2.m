% Calculate vector of cutting forces with Newton-Euler
% S5PRPRR8
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-12-05 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRPRR8_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR8_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR8_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR8_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR8_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR8_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR8_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:02:31
% EndTime: 2019-12-05 16:02:32
% DurationCPUTime: 0.57s
% Computational Cost: add. (5213->109), mult. (9445->145), div. (0->0), fcn. (5830->10), ass. (0->68)
t56 = sin(pkin(9));
t58 = cos(pkin(9));
t47 = -t58 * g(1) - t56 * g(2);
t55 = -g(3) + qJDD(1);
t62 = sin(qJ(2));
t65 = cos(qJ(2));
t57 = sin(pkin(5));
t89 = t57 * t62;
t46 = t56 * g(1) - t58 * g(2);
t59 = cos(pkin(5));
t90 = t46 * t59;
t80 = t65 * t47 + t55 * t89 + t62 * t90;
t94 = -qJDD(2) * qJ(3) - (2 * qJD(3) * qJD(2)) - t80;
t93 = -t62 * t47 + (t55 * t57 + t90) * t65;
t92 = -pkin(2) - pkin(7);
t67 = qJD(2) ^ 2;
t61 = sin(qJ(4));
t64 = cos(qJ(4));
t43 = (pkin(4) * t61 - pkin(8) * t64) * qJD(2);
t66 = qJD(4) ^ 2;
t83 = t61 * qJD(2);
t69 = -t67 * qJ(3) + qJDD(3) - t93;
t22 = t92 * qJDD(2) + t69;
t32 = -t57 * t46 + t59 * t55;
t85 = t61 * t22 + t64 * t32;
t17 = -t66 * pkin(4) + qJDD(4) * pkin(8) - t43 * t83 + t85;
t82 = qJD(2) * qJD(4);
t77 = t64 * t82;
t44 = -t61 * qJDD(2) - t77;
t78 = t61 * t82;
t45 = t64 * qJDD(2) - t78;
t71 = t92 * t67 - t94;
t18 = (-t45 + t78) * pkin(8) + (-t44 + t77) * pkin(4) + t71;
t60 = sin(qJ(5));
t63 = cos(qJ(5));
t84 = qJD(2) * t64;
t40 = t63 * qJD(4) - t60 * t84;
t26 = t40 * qJD(5) + t60 * qJDD(4) + t63 * t45;
t41 = t60 * qJD(4) + t63 * t84;
t27 = -t40 * mrSges(6,1) + t41 * mrSges(6,2);
t51 = qJD(5) + t83;
t29 = -t51 * mrSges(6,2) + t40 * mrSges(6,3);
t37 = qJDD(5) - t44;
t14 = m(6) * (-t60 * t17 + t63 * t18) - t26 * mrSges(6,3) + t37 * mrSges(6,1) - t41 * t27 + t51 * t29;
t25 = -t41 * qJD(5) + t63 * qJDD(4) - t60 * t45;
t30 = t51 * mrSges(6,1) - t41 * mrSges(6,3);
t15 = m(6) * (t63 * t17 + t60 * t18) + t25 * mrSges(6,3) - t37 * mrSges(6,2) + t40 * t27 - t51 * t30;
t42 = (mrSges(5,1) * t61 + mrSges(5,2) * t64) * qJD(2);
t49 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t84;
t10 = m(5) * t85 - qJDD(4) * mrSges(5,2) + t44 * mrSges(5,3) - qJD(4) * t49 - t60 * t14 + t63 * t15 - t42 * t83;
t48 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t83;
t88 = t61 * t32;
t68 = m(6) * (-qJDD(4) * pkin(4) - t66 * pkin(8) + t88 + (qJD(2) * t43 - t22) * t64) - t25 * mrSges(6,1) + t26 * mrSges(6,2) - t40 * t29 + t41 * t30;
t11 = m(5) * (t64 * t22 - t88) - t45 * mrSges(5,3) + qJDD(4) * mrSges(5,1) - t42 * t84 + qJD(4) * t48 - t68;
t73 = -m(4) * (-qJDD(2) * pkin(2) + t69) - t61 * t10 - t64 * t11;
t86 = (-mrSges(3,2) + mrSges(4,3));
t87 = mrSges(3,1) - mrSges(4,2);
t4 = m(3) * t93 + t87 * qJDD(2) + (t86 * t67) + t73;
t91 = t4 * t65;
t75 = m(4) * t32 + t64 * t10 - t61 * t11;
t6 = m(3) * t32 + t75;
t72 = m(5) * t71 - t44 * mrSges(5,1) + t45 * mrSges(5,2) + t63 * t14 + t60 * t15 + t48 * t83 + t49 * t84;
t70 = -m(4) * (t67 * pkin(2) + t94) + t72;
t8 = m(3) * t80 + t86 * qJDD(2) - t87 * t67 + t70;
t79 = m(2) * t55 + t57 * t91 + t59 * t6 + t8 * t89;
t2 = m(2) * t47 - t62 * t4 + t65 * t8;
t1 = m(2) * t46 - t57 * t6 + (t62 * t8 + t91) * t59;
t3 = [-m(1) * g(1) - t56 * t1 + t58 * t2, t2, t8, t75, t10, t15; -m(1) * g(2) + t58 * t1 + t56 * t2, t1, t4, -(t67 * mrSges(4,2)) - qJDD(2) * mrSges(4,3) - t70, t11, t14; -m(1) * g(3) + t79, t79, t6, qJDD(2) * mrSges(4,2) - t67 * mrSges(4,3) - t73, t72, t68;];
f_new = t3;
