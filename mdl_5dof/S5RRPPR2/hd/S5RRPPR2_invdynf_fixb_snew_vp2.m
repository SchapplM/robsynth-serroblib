% Calculate vector of cutting forces with Newton-Euler
% S5RRPPR2
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
% Datum: 2022-01-20 10:06
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRPPR2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR2_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR2_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR2_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:05:33
% EndTime: 2022-01-20 10:05:35
% DurationCPUTime: 0.70s
% Computational Cost: add. (8897->100), mult. (11760->136), div. (0->0), fcn. (6518->10), ass. (0->65)
t48 = sin(pkin(9));
t50 = cos(pkin(9));
t83 = (t48 ^ 2 + t50 ^ 2) * mrSges(5,3);
t82 = 2 * qJD(4);
t81 = -m(2) - m(3);
t43 = qJDD(1) + qJDD(2);
t80 = t43 * mrSges(5,3);
t46 = qJD(1) + qJD(2);
t79 = t46 * t48;
t78 = t50 * t46;
t47 = -g(3) + qJDD(3);
t77 = t50 * t47;
t54 = sin(qJ(1));
t57 = cos(qJ(1));
t68 = t54 * g(1) - t57 * g(2);
t34 = qJDD(1) * pkin(1) + t68;
t58 = qJD(1) ^ 2;
t64 = -t57 * g(1) - t54 * g(2);
t35 = -t58 * pkin(1) + t64;
t53 = sin(qJ(2));
t56 = cos(qJ(2));
t66 = t56 * t34 - t53 * t35;
t22 = t43 * pkin(2) + t66;
t42 = t46 ^ 2;
t75 = t53 * t34 + t56 * t35;
t23 = -t42 * pkin(2) + t75;
t49 = sin(pkin(8));
t51 = cos(pkin(8));
t76 = t49 * t22 + t51 * t23;
t73 = qJD(5) * t46;
t18 = -t42 * pkin(3) + t43 * qJ(4) + t76;
t63 = -t50 * mrSges(5,1) + t48 * mrSges(5,2);
t29 = t63 * t46;
t52 = sin(qJ(5));
t55 = cos(qJ(5));
t27 = (-t43 * t52 - t55 * t73) * t48;
t28 = (t43 * t55 - t52 * t73) * t48;
t65 = -pkin(4) * t50 - pkin(7) * t48;
t30 = t65 * t46;
t60 = m(6) * (-t77 + (t18 + (t82 + t30) * t46) * t48) - t27 * mrSges(6,1) + t28 * mrSges(6,2);
t37 = qJD(5) - t78;
t71 = mrSges(6,3) * t79;
t24 = -t37 * mrSges(6,2) - t52 * t71;
t25 = t37 * mrSges(6,1) - t55 * t71;
t62 = t52 * t24 + t55 * t25;
t10 = m(5) * t77 + (-m(5) * t18 - t80 + (-0.2e1 * m(5) * qJD(4) - t29 - t62) * t46) * t48 - t60;
t69 = t50 * t18 + t48 * t47 + t78 * t82;
t14 = t30 * t78 + t69;
t67 = t51 * t22 - t49 * t23;
t59 = -t42 * qJ(4) + qJDD(4) - t67;
t15 = (-pkin(3) + t65) * t43 + t59;
t36 = -t50 * t43 + qJDD(5);
t70 = (mrSges(6,1) * t52 + mrSges(6,2) * t55) * t79 ^ 2;
t11 = m(6) * (-t52 * t14 + t55 * t15) - t28 * mrSges(6,3) + t36 * mrSges(6,1) - t55 * t70 + t37 * t24;
t12 = m(6) * (t55 * t14 + t52 * t15) + t27 * mrSges(6,3) - t36 * mrSges(6,2) - t52 * t70 - t37 * t25;
t8 = m(5) * t69 + t55 * t12 - t52 * t11 + (t46 * t29 + t80) * t50;
t72 = m(4) * t47 + t50 * t10 + t48 * t8;
t61 = m(5) * (-t43 * pkin(3) + t59) + t55 * t11 + t52 * t12;
t6 = m(4) * t67 + (mrSges(4,1) - t63) * t43 + (-mrSges(4,2) + t83) * t42 - t61;
t5 = m(4) * t76 - t42 * mrSges(4,1) - t43 * mrSges(4,2) - t48 * t10 + t50 * t8;
t4 = m(3) * t75 - t42 * mrSges(3,1) - t43 * mrSges(3,2) - t49 * t6 + t51 * t5;
t3 = m(3) * t66 + t43 * mrSges(3,1) - t42 * mrSges(3,2) + t49 * t5 + t51 * t6;
t2 = m(2) * t64 - t58 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t53 * t3 + t56 * t4;
t1 = m(2) * t68 + qJDD(1) * mrSges(2,1) - t58 * mrSges(2,2) + t56 * t3 + t53 * t4;
t7 = [-m(1) * g(1) - t54 * t1 + t57 * t2, t2, t4, t5, t8, t12; -m(1) * g(2) + t57 * t1 + t54 * t2, t1, t3, t6, t10, t11; (-m(1) + t81) * g(3) + t72, t81 * g(3) + t72, -m(3) * g(3) + t72, t72, -t42 * t83 + t63 * t43 + t61, t62 * t79 + t60;];
f_new = t7;
