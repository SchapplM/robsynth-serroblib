% Calculate vector of cutting forces with Newton-Euler
% S5RRRPR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2022-01-20 11:31
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRPR2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR2_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR2_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR2_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR2_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:30:38
% EndTime: 2022-01-20 11:30:38
% DurationCPUTime: 0.67s
% Computational Cost: add. (10542->86), mult. (12295->114), div. (0->0), fcn. (6580->10), ass. (0->56)
t40 = qJDD(1) + qJDD(2);
t35 = qJDD(3) + t40;
t45 = sin(qJ(5));
t49 = cos(qJ(5));
t41 = qJD(1) + qJD(2);
t36 = qJD(3) + t41;
t63 = qJD(5) * t36;
t25 = t45 * t35 + t49 * t63;
t26 = t49 * t35 - t45 * t63;
t68 = t36 * t45;
t30 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t68;
t67 = t36 * t49;
t31 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t67;
t34 = t36 ^ 2;
t48 = sin(qJ(1));
t52 = cos(qJ(1));
t60 = t48 * g(1) - t52 * g(2);
t32 = qJDD(1) * pkin(1) + t60;
t53 = qJD(1) ^ 2;
t57 = -t52 * g(1) - t48 * g(2);
t33 = -t53 * pkin(1) + t57;
t47 = sin(qJ(2));
t51 = cos(qJ(2));
t58 = t51 * t32 - t47 * t33;
t22 = t40 * pkin(2) + t58;
t39 = t41 ^ 2;
t64 = t47 * t32 + t51 * t33;
t23 = -t39 * pkin(2) + t64;
t46 = sin(qJ(3));
t50 = cos(qJ(3));
t59 = t50 * t22 - t46 * t23;
t17 = t35 * pkin(3) + t59;
t65 = t46 * t22 + t50 * t23;
t18 = -t34 * pkin(3) + t65;
t43 = sin(pkin(9));
t44 = cos(pkin(9));
t56 = t44 * t17 - t43 * t18;
t70 = (t45 * t30 - t49 * t31) * t36 + m(6) * (-t35 * pkin(4) - t34 * pkin(8) - t56) - t26 * mrSges(6,1) + t25 * mrSges(6,2);
t69 = -m(3) - m(4);
t66 = t43 * t17 + t44 * t18;
t62 = -m(2) + t69;
t14 = -t34 * pkin(4) + t35 * pkin(8) + t66;
t24 = (-mrSges(6,1) * t49 + mrSges(6,2) * t45) * t36;
t42 = -g(3) + qJDD(4);
t11 = m(6) * (-t45 * t14 + t49 * t42) - t25 * mrSges(6,3) + qJDD(5) * mrSges(6,1) - t24 * t68 + qJD(5) * t31;
t12 = m(6) * (t49 * t14 + t45 * t42) + t26 * mrSges(6,3) - qJDD(5) * mrSges(6,2) + t24 * t67 - qJD(5) * t30;
t61 = m(5) * t42 + t49 * t11 + t45 * t12;
t8 = m(5) * t56 + t35 * mrSges(5,1) - t34 * mrSges(5,2) - t70;
t7 = m(5) * t66 - t34 * mrSges(5,1) - t35 * mrSges(5,2) - t45 * t11 + t49 * t12;
t6 = m(4) * t65 - t34 * mrSges(4,1) - t35 * mrSges(4,2) - t43 * t8 + t44 * t7;
t5 = m(4) * t59 + t35 * mrSges(4,1) - t34 * mrSges(4,2) + t43 * t7 + t44 * t8;
t4 = m(3) * t64 - t39 * mrSges(3,1) - t40 * mrSges(3,2) - t46 * t5 + t50 * t6;
t3 = m(3) * t58 + t40 * mrSges(3,1) - t39 * mrSges(3,2) + t46 * t6 + t50 * t5;
t2 = m(2) * t57 - t53 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t47 * t3 + t51 * t4;
t1 = m(2) * t60 + qJDD(1) * mrSges(2,1) - t53 * mrSges(2,2) + t51 * t3 + t47 * t4;
t9 = [-m(1) * g(1) - t48 * t1 + t52 * t2, t2, t4, t6, t7, t12; -m(1) * g(2) + t52 * t1 + t48 * t2, t1, t3, t5, t8, t11; (-m(1) + t62) * g(3) + t61, t62 * g(3) + t61, t69 * g(3) + t61, -m(4) * g(3) + t61, t61, t70;];
f_new = t9;
