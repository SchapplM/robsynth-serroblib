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
% Datum: 2019-12-05 18:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:40:49
% EndTime: 2019-12-05 18:40:49
% DurationCPUTime: 0.60s
% Computational Cost: add. (10542->86), mult. (12295->114), div. (0->0), fcn. (6580->10), ass. (0->56)
t42 = qJDD(1) + qJDD(2);
t35 = qJDD(3) + t42;
t47 = sin(qJ(5));
t51 = cos(qJ(5));
t43 = qJD(1) + qJD(2);
t36 = qJD(3) + t43;
t64 = qJD(5) * t36;
t25 = t47 * t35 + t51 * t64;
t26 = t51 * t35 - t47 * t64;
t70 = t36 * t47;
t30 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t70;
t69 = t36 * t51;
t31 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t69;
t34 = t36 ^ 2;
t50 = sin(qJ(1));
t54 = cos(qJ(1));
t65 = t54 * g(2) + t50 * g(3);
t32 = qJDD(1) * pkin(1) + t65;
t55 = qJD(1) ^ 2;
t61 = t50 * g(2) - t54 * g(3);
t33 = -t55 * pkin(1) + t61;
t49 = sin(qJ(2));
t53 = cos(qJ(2));
t59 = t53 * t32 - t49 * t33;
t22 = t42 * pkin(2) + t59;
t41 = t43 ^ 2;
t66 = t49 * t32 + t53 * t33;
t23 = -t41 * pkin(2) + t66;
t48 = sin(qJ(3));
t52 = cos(qJ(3));
t60 = t52 * t22 - t48 * t23;
t17 = t35 * pkin(3) + t60;
t67 = t48 * t22 + t52 * t23;
t18 = -t34 * pkin(3) + t67;
t45 = sin(pkin(9));
t46 = cos(pkin(9));
t58 = t46 * t17 - t45 * t18;
t72 = (t47 * t30 - t51 * t31) * t36 + m(6) * (-t35 * pkin(4) - t34 * pkin(8) - t58) - t26 * mrSges(6,1) + t25 * mrSges(6,2);
t71 = -m(3) - m(4);
t68 = t45 * t17 + t46 * t18;
t63 = -m(2) + t71;
t14 = -t34 * pkin(4) + t35 * pkin(8) + t68;
t24 = (-mrSges(6,1) * t51 + mrSges(6,2) * t47) * t36;
t44 = -g(1) + qJDD(4);
t11 = m(6) * (-t47 * t14 + t51 * t44) - t25 * mrSges(6,3) + qJDD(5) * mrSges(6,1) - t24 * t70 + qJD(5) * t31;
t12 = m(6) * (t51 * t14 + t47 * t44) + t26 * mrSges(6,3) - qJDD(5) * mrSges(6,2) + t24 * t69 - qJD(5) * t30;
t62 = m(5) * t44 + t51 * t11 + t47 * t12;
t8 = m(5) * t58 + t35 * mrSges(5,1) - t34 * mrSges(5,2) - t72;
t7 = m(5) * t68 - t34 * mrSges(5,1) - t35 * mrSges(5,2) - t47 * t11 + t51 * t12;
t6 = m(4) * t67 - t34 * mrSges(4,1) - t35 * mrSges(4,2) - t45 * t8 + t46 * t7;
t5 = m(4) * t60 + t35 * mrSges(4,1) - t34 * mrSges(4,2) + t45 * t7 + t46 * t8;
t4 = m(3) * t66 - t41 * mrSges(3,1) - t42 * mrSges(3,2) - t48 * t5 + t52 * t6;
t3 = m(3) * t59 + t42 * mrSges(3,1) - t41 * mrSges(3,2) + t48 * t6 + t52 * t5;
t2 = m(2) * t65 + qJDD(1) * mrSges(2,1) - t55 * mrSges(2,2) + t53 * t3 + t49 * t4;
t1 = m(2) * t61 - t55 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t49 * t3 + t53 * t4;
t9 = [(-m(1) + t63) * g(1) + t62, t1, t4, t6, t7, t12; -m(1) * g(2) - t50 * t1 - t54 * t2, t2, t3, t5, t8, t11; -m(1) * g(3) + t54 * t1 - t50 * t2, t63 * g(1) + t62, t71 * g(1) + t62, -m(4) * g(1) + t62, t62, t72;];
f_new = t9;
