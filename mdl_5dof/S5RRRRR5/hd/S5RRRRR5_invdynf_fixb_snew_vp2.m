% Calculate vector of cutting forces with Newton-Euler
% S5RRRRR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2022-01-20 12:02
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRRR5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR5_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR5_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR5_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR5_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 12:01:47
% EndTime: 2022-01-20 12:01:50
% DurationCPUTime: 0.89s
% Computational Cost: add. (14593->116), mult. (15052->154), div. (0->0), fcn. (8503->10), ass. (0->69)
t51 = qJDD(1) + qJDD(2);
t46 = qJDD(3) + t51;
t56 = sin(qJ(4));
t61 = cos(qJ(4));
t53 = qJD(1) + qJD(2);
t47 = qJD(3) + t53;
t76 = qJD(4) * t47;
t35 = t56 * t46 + t61 * t76;
t36 = t61 * t46 - t56 * t76;
t81 = t47 * t56;
t40 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t81;
t80 = t47 * t61;
t41 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t80;
t45 = t47 ^ 2;
t55 = sin(qJ(5));
t60 = cos(qJ(5));
t33 = (t55 * t61 + t56 * t60) * t47;
t21 = -t33 * qJD(5) - t55 * t35 + t60 * t36;
t32 = (-t55 * t56 + t60 * t61) * t47;
t22 = t32 * qJD(5) + t60 * t35 + t55 * t36;
t52 = qJD(4) + qJD(5);
t30 = -t52 * mrSges(6,2) + t32 * mrSges(6,3);
t31 = t52 * mrSges(6,1) - t33 * mrSges(6,3);
t42 = qJD(4) * pkin(4) - pkin(9) * t81;
t54 = t61 ^ 2;
t59 = sin(qJ(1));
t64 = cos(qJ(1));
t74 = t59 * g(1) - t64 * g(2);
t43 = qJDD(1) * pkin(1) + t74;
t65 = qJD(1) ^ 2;
t70 = -t64 * g(1) - t59 * g(2);
t44 = -t65 * pkin(1) + t70;
t58 = sin(qJ(2));
t63 = cos(qJ(2));
t71 = t63 * t43 - t58 * t44;
t28 = t51 * pkin(2) + t71;
t49 = t53 ^ 2;
t77 = t58 * t43 + t63 * t44;
t29 = -t49 * pkin(2) + t77;
t57 = sin(qJ(3));
t62 = cos(qJ(3));
t72 = t62 * t28 - t57 * t29;
t68 = -t46 * pkin(3) - t72;
t67 = t21 * mrSges(6,1) + t32 * t30 - m(6) * (t42 * t81 - t36 * pkin(4) + (-pkin(9) * t54 - pkin(8)) * t45 + t68) - t22 * mrSges(6,2) - t33 * t31;
t85 = (t56 * t40 - t61 * t41) * t47 + m(5) * (-t45 * pkin(8) + t68) - t36 * mrSges(5,1) + t35 * mrSges(5,2) - t67;
t84 = -m(3) - m(4);
t78 = t57 * t28 + t62 * t29;
t19 = -t45 * pkin(3) + t46 * pkin(8) + t78;
t79 = t56 * t19;
t82 = pkin(4) * t45;
t13 = qJDD(4) * pkin(4) - t35 * pkin(9) - t79 + (pkin(9) * t76 + t56 * t82 - g(3)) * t61;
t73 = -t56 * g(3) + t61 * t19;
t14 = t36 * pkin(9) - qJD(4) * t42 - t54 * t82 + t73;
t24 = -t32 * mrSges(6,1) + t33 * mrSges(6,2);
t50 = qJDD(4) + qJDD(5);
t11 = m(6) * (t60 * t13 - t55 * t14) - t22 * mrSges(6,3) + t50 * mrSges(6,1) - t33 * t24 + t52 * t30;
t12 = m(6) * (t55 * t13 + t60 * t14) + t21 * mrSges(6,3) - t50 * mrSges(6,2) + t32 * t24 - t52 * t31;
t34 = (-mrSges(5,1) * t61 + mrSges(5,2) * t56) * t47;
t8 = m(5) * (-t61 * g(3) - t79) - t35 * mrSges(5,3) + qJDD(4) * mrSges(5,1) - t34 * t81 + qJD(4) * t41 + t55 * t12 + t60 * t11;
t9 = m(5) * t73 - qJDD(4) * mrSges(5,2) + t36 * mrSges(5,3) - qJD(4) * t40 - t55 * t11 + t60 * t12 + t34 * t80;
t83 = t56 * t9 + t61 * t8;
t75 = -m(2) + t84;
t10 = m(4) * t72 + t46 * mrSges(4,1) - t45 * mrSges(4,2) - t85;
t5 = m(4) * t78 - t45 * mrSges(4,1) - t46 * mrSges(4,2) - t56 * t8 + t61 * t9;
t4 = m(3) * t77 - t49 * mrSges(3,1) - t51 * mrSges(3,2) - t57 * t10 + t62 * t5;
t3 = m(3) * t71 + t51 * mrSges(3,1) - t49 * mrSges(3,2) + t62 * t10 + t57 * t5;
t2 = m(2) * t70 - t65 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t58 * t3 + t63 * t4;
t1 = m(2) * t74 + qJDD(1) * mrSges(2,1) - t65 * mrSges(2,2) + t63 * t3 + t58 * t4;
t6 = [-m(1) * g(1) - t59 * t1 + t64 * t2, t2, t4, t5, t9, t12; -m(1) * g(2) + t64 * t1 + t59 * t2, t1, t3, t10, t8, t11; (-m(1) + t75) * g(3) + t83, t75 * g(3) + t83, t84 * g(3) + t83, -m(4) * g(3) + t83, t85, -t67;];
f_new = t6;
