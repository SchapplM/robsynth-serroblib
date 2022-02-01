% Calculate vector of cutting forces with Newton-Euler
% S5RRPPR1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2022-01-20 09:52
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRPPR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR1_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR1_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR1_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR1_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:51:25
% EndTime: 2022-01-20 09:51:26
% DurationCPUTime: 0.75s
% Computational Cost: add. (10111->101), mult. (13789->133), div. (0->0), fcn. (7969->10), ass. (0->63)
t50 = qJD(1) + qJD(2);
t46 = t50 ^ 2;
t54 = cos(pkin(9));
t49 = t54 ^ 2;
t52 = sin(pkin(9));
t77 = t52 ^ 2 + t49;
t84 = t77 * mrSges(5,3);
t83 = -m(2) - m(3);
t82 = pkin(4) * t54;
t47 = qJDD(1) + qJDD(2);
t81 = pkin(7) * t47;
t58 = sin(qJ(1));
t61 = cos(qJ(1));
t73 = t58 * g(1) - t61 * g(2);
t38 = qJDD(1) * pkin(1) + t73;
t62 = qJD(1) ^ 2;
t70 = -t61 * g(1) - t58 * g(2);
t39 = -t62 * pkin(1) + t70;
t57 = sin(qJ(2));
t60 = cos(qJ(2));
t71 = t60 * t38 - t57 * t39;
t28 = t47 * pkin(2) + t71;
t79 = t57 * t38 + t60 * t39;
t29 = -t46 * pkin(2) + t79;
t53 = sin(pkin(8));
t55 = cos(pkin(8));
t80 = t53 * t28 + t55 * t29;
t51 = -g(3) + qJDD(3);
t76 = qJD(4) * t50;
t78 = t54 * t51 - 0.2e1 * t52 * t76;
t19 = -t46 * pkin(3) + t47 * qJ(4) + t80;
t13 = (t46 * t82 - t19 - t81) * t52 + t78;
t74 = t52 * t51 + (t19 + 0.2e1 * t76) * t54;
t14 = -t49 * t46 * pkin(4) + t54 * t81 + t74;
t56 = sin(qJ(5));
t59 = cos(qJ(5));
t65 = -t52 * t56 + t54 * t59;
t32 = t65 * t50;
t66 = t52 * t59 + t54 * t56;
t33 = t66 * t50;
t22 = -t32 * mrSges(6,1) + t33 * mrSges(6,2);
t24 = t32 * qJD(5) + t66 * t47;
t30 = -qJD(5) * mrSges(6,2) + t32 * mrSges(6,3);
t11 = m(6) * (t59 * t13 - t56 * t14) - t24 * mrSges(6,3) + qJDD(5) * mrSges(6,1) - t33 * t22 + qJD(5) * t30;
t23 = -t33 * qJD(5) + t65 * t47;
t31 = qJD(5) * mrSges(6,1) - t33 * mrSges(6,3);
t12 = m(6) * (t56 * t13 + t59 * t14) + t23 * mrSges(6,3) - qJDD(5) * mrSges(6,2) + t32 * t22 - qJD(5) * t31;
t68 = -t54 * mrSges(5,1) + t52 * mrSges(5,2);
t67 = t47 * mrSges(5,3) + t46 * t68;
t8 = m(5) * t78 + t56 * t12 + t59 * t11 + (-m(5) * t19 - t67) * t52;
t9 = m(5) * t74 - t56 * t11 + t59 * t12 + t67 * t54;
t75 = m(4) * t51 + t52 * t9 + t54 * t8;
t72 = t55 * t28 - t53 * t29;
t69 = qJDD(4) - t72;
t64 = t23 * mrSges(6,1) + t32 * t30 - m(6) * ((-pkin(3) - t82) * t47 + (-t77 * pkin(7) - qJ(4)) * t46 + t69) - t33 * t31 - t24 * mrSges(6,2);
t63 = m(5) * (-t47 * pkin(3) - t46 * qJ(4) + t69) - t64;
t10 = m(4) * t72 + (mrSges(4,1) - t68) * t47 + (-mrSges(4,2) + t84) * t46 - t63;
t5 = m(4) * t80 - t46 * mrSges(4,1) - t47 * mrSges(4,2) - t52 * t8 + t54 * t9;
t4 = m(3) * t79 - t46 * mrSges(3,1) - t47 * mrSges(3,2) - t53 * t10 + t55 * t5;
t3 = m(3) * t71 + t47 * mrSges(3,1) - t46 * mrSges(3,2) + t55 * t10 + t53 * t5;
t2 = m(2) * t70 - t62 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t57 * t3 + t60 * t4;
t1 = m(2) * t73 + qJDD(1) * mrSges(2,1) - t62 * mrSges(2,2) + t60 * t3 + t57 * t4;
t6 = [-m(1) * g(1) - t58 * t1 + t61 * t2, t2, t4, t5, t9, t12; -m(1) * g(2) + t61 * t1 + t58 * t2, t1, t3, t10, t8, t11; (-m(1) + t83) * g(3) + t75, t83 * g(3) + t75, -m(3) * g(3) + t75, t75, -t46 * t84 + t68 * t47 + t63, -t64;];
f_new = t6;
