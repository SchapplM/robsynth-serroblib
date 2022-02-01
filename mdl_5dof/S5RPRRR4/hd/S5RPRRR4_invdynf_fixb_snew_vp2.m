% Calculate vector of cutting forces with Newton-Euler
% S5RPRRR4
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
% Datum: 2022-01-23 09:35
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRRR4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR4_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR4_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR4_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR4_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:34:32
% EndTime: 2022-01-23 09:34:33
% DurationCPUTime: 0.61s
% Computational Cost: add. (9393->84), mult. (12295->114), div. (0->0), fcn. (6580->10), ass. (0->56)
t42 = qJDD(1) + qJDD(3);
t35 = qJDD(4) + t42;
t47 = sin(qJ(5));
t51 = cos(qJ(5));
t43 = qJD(1) + qJD(3);
t36 = qJD(4) + t43;
t66 = qJD(5) * t36;
t25 = t47 * t35 + t51 * t66;
t26 = t51 * t35 - t47 * t66;
t71 = t36 * t47;
t30 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t71;
t70 = t36 * t51;
t31 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t70;
t34 = t36 ^ 2;
t50 = sin(qJ(1));
t54 = cos(qJ(1));
t64 = t50 * g(1) - t54 * g(2);
t32 = qJDD(1) * pkin(1) + t64;
t55 = qJD(1) ^ 2;
t59 = -t54 * g(1) - t50 * g(2);
t33 = -t55 * pkin(1) + t59;
t45 = sin(pkin(9));
t46 = cos(pkin(9));
t62 = t46 * t32 - t45 * t33;
t22 = qJDD(1) * pkin(2) + t62;
t67 = t45 * t32 + t46 * t33;
t23 = -t55 * pkin(2) + t67;
t49 = sin(qJ(3));
t53 = cos(qJ(3));
t63 = t53 * t22 - t49 * t23;
t17 = t42 * pkin(3) + t63;
t41 = t43 ^ 2;
t68 = t49 * t22 + t53 * t23;
t18 = -t41 * pkin(3) + t68;
t48 = sin(qJ(4));
t52 = cos(qJ(4));
t58 = t52 * t17 - t48 * t18;
t72 = (t47 * t30 - t51 * t31) * t36 + m(6) * (-t35 * pkin(4) - t34 * pkin(8) - t58) - t26 * mrSges(6,1) + t25 * mrSges(6,2);
t69 = t48 * t17 + t52 * t18;
t14 = -t34 * pkin(4) + t35 * pkin(8) + t69;
t24 = (-mrSges(6,1) * t51 + mrSges(6,2) * t47) * t36;
t44 = -g(3) + qJDD(2);
t11 = m(6) * (-t47 * t14 + t51 * t44) - t25 * mrSges(6,3) + qJDD(5) * mrSges(6,1) - t24 * t71 + qJD(5) * t31;
t12 = m(6) * (t51 * t14 + t47 * t44) + t26 * mrSges(6,3) - qJDD(5) * mrSges(6,2) + t24 * t70 - qJD(5) * t30;
t65 = m(5) * t44 + t51 * t11 + t47 * t12;
t61 = m(4) * t44 + t65;
t60 = m(3) * t44 + t61;
t8 = m(5) * t58 + t35 * mrSges(5,1) - t34 * mrSges(5,2) - t72;
t7 = m(5) * t69 - t34 * mrSges(5,1) - t35 * mrSges(5,2) - t47 * t11 + t51 * t12;
t6 = m(4) * t68 - t41 * mrSges(4,1) - t42 * mrSges(4,2) - t48 * t8 + t52 * t7;
t5 = m(4) * t63 + t42 * mrSges(4,1) - t41 * mrSges(4,2) + t48 * t7 + t52 * t8;
t4 = m(3) * t67 - t55 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t49 * t5 + t53 * t6;
t3 = m(3) * t62 + qJDD(1) * mrSges(3,1) - t55 * mrSges(3,2) + t49 * t6 + t53 * t5;
t2 = m(2) * t59 - t55 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t45 * t3 + t46 * t4;
t1 = m(2) * t64 + qJDD(1) * mrSges(2,1) - t55 * mrSges(2,2) + t46 * t3 + t45 * t4;
t9 = [-m(1) * g(1) - t50 * t1 + t54 * t2, t2, t4, t6, t7, t12; -m(1) * g(2) + t54 * t1 + t50 * t2, t1, t3, t5, t8, t11; (-m(1) - m(2)) * g(3) + t60, -m(2) * g(3) + t60, t60, t61, t65, t72;];
f_new = t9;
