% Calculate vector of cutting forces with Newton-Euler
% S5PPRRR2
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
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PPRRR2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR2_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR2_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR2_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR2_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR2_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR2_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR2_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:14:15
% EndTime: 2019-12-05 15:14:16
% DurationCPUTime: 0.60s
% Computational Cost: add. (6115->98), mult. (11125->139), div. (0->0), fcn. (7335->10), ass. (0->61)
t57 = sin(qJ(4));
t60 = cos(qJ(4));
t74 = qJD(3) * qJD(4);
t71 = t60 * t74;
t40 = qJDD(3) * t57 + t71;
t41 = qJDD(3) * t60 - t57 * t74;
t76 = qJD(3) * t57;
t44 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t76;
t75 = qJD(3) * t60;
t45 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t75;
t62 = qJD(3) ^ 2;
t56 = sin(qJ(5));
t59 = cos(qJ(5));
t35 = (t56 * t60 + t57 * t59) * qJD(3);
t23 = -qJD(5) * t35 - t40 * t56 + t41 * t59;
t34 = (-t56 * t57 + t59 * t60) * qJD(3);
t24 = qJD(5) * t34 + t40 * t59 + t41 * t56;
t49 = qJD(4) + qJD(5);
t30 = -mrSges(6,2) * t49 + mrSges(6,3) * t34;
t31 = mrSges(6,1) * t49 - mrSges(6,3) * t35;
t46 = qJD(4) * pkin(4) - pkin(7) * t76;
t50 = t60 ^ 2;
t53 = sin(pkin(8));
t55 = cos(pkin(8));
t43 = -g(1) * t55 - g(2) * t53;
t51 = -g(3) + qJDD(1);
t52 = sin(pkin(9));
t54 = cos(pkin(9));
t32 = -t43 * t52 + t51 * t54;
t33 = t43 * t54 + t51 * t52;
t58 = sin(qJ(3));
t61 = cos(qJ(3));
t69 = t32 * t61 - t58 * t33;
t65 = -qJDD(3) * pkin(3) - t69;
t64 = t23 * mrSges(6,1) + t34 * t30 - m(6) * (t46 * t76 - pkin(4) * t41 + (-pkin(7) * t50 - pkin(6)) * t62 + t65) - t24 * mrSges(6,2) - t35 * t31;
t79 = (t44 * t57 - t45 * t60) * qJD(3) + m(5) * (-pkin(6) * t62 + t65) - t41 * mrSges(5,1) + t40 * mrSges(5,2) - t64;
t77 = t58 * t32 + t61 * t33;
t21 = -pkin(3) * t62 + qJDD(3) * pkin(6) + t77;
t68 = g(1) * t53 - g(2) * t55;
t42 = qJDD(2) - t68;
t78 = t60 * t21 + t57 * t42;
t70 = -t57 * t21 + t60 * t42;
t15 = (-t40 + t71) * pkin(7) + (t57 * t60 * t62 + qJDD(4)) * pkin(4) + t70;
t16 = -pkin(4) * t50 * t62 + pkin(7) * t41 - qJD(4) * t46 + t78;
t26 = -mrSges(6,1) * t34 + mrSges(6,2) * t35;
t48 = qJDD(4) + qJDD(5);
t13 = m(6) * (t15 * t59 - t16 * t56) - t24 * mrSges(6,3) + t48 * mrSges(6,1) - t35 * t26 + t49 * t30;
t14 = m(6) * (t15 * t56 + t16 * t59) + t23 * mrSges(6,3) - t48 * mrSges(6,2) + t34 * t26 - t49 * t31;
t39 = (-mrSges(5,1) * t60 + mrSges(5,2) * t57) * qJD(3);
t10 = m(5) * t70 + qJDD(4) * mrSges(5,1) - t40 * mrSges(5,3) + qJD(4) * t45 + t59 * t13 + t56 * t14 - t39 * t76;
t11 = m(5) * t78 - qJDD(4) * mrSges(5,2) + t41 * mrSges(5,3) - qJD(4) * t44 - t56 * t13 + t59 * t14 + t39 * t75;
t73 = m(4) * t42 + t60 * t10 + t57 * t11;
t12 = m(4) * t69 + qJDD(3) * mrSges(4,1) - t62 * mrSges(4,2) - t79;
t6 = m(4) * t77 - t62 * mrSges(4,1) - qJDD(3) * mrSges(4,2) - t57 * t10 + t60 * t11;
t4 = m(3) * t32 + t12 * t61 + t58 * t6;
t5 = m(3) * t33 - t12 * t58 + t6 * t61;
t72 = m(2) * t51 + t54 * t4 + t52 * t5;
t67 = m(3) * t42 + t73;
t7 = m(2) * t68 - t67;
t1 = m(2) * t43 - t4 * t52 + t5 * t54;
t2 = [-m(1) * g(1) + t1 * t55 - t53 * t7, t1, t5, t6, t11, t14; -m(1) * g(2) + t1 * t53 + t55 * t7, t7, t4, t12, t10, t13; -m(1) * g(3) + t72, t72, t67, t73, t79, -t64;];
f_new = t2;
