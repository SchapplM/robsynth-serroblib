% Calculate vector of cutting forces with Newton-Euler
% S5PRRPR8
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRRPR8_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR8_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR8_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR8_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR8_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR8_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR8_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:42:22
% EndTime: 2019-12-31 17:42:23
% DurationCPUTime: 0.47s
% Computational Cost: add. (6209->78), mult. (7808->107), div. (0->0), fcn. (4710->10), ass. (0->52)
t42 = qJDD(2) + qJDD(3);
t48 = sin(qJ(5));
t51 = cos(qJ(5));
t43 = qJD(2) + qJD(3);
t63 = qJD(5) * t43;
t27 = t48 * t42 + t51 * t63;
t28 = t51 * t42 - t48 * t63;
t69 = t43 * t48;
t31 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t69;
t68 = t43 * t51;
t32 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t68;
t41 = t43 ^ 2;
t46 = sin(pkin(8));
t64 = cos(pkin(8));
t36 = -t64 * g(1) - t46 * g(2);
t44 = -g(3) + qJDD(1);
t50 = sin(qJ(2));
t53 = cos(qJ(2));
t59 = -t50 * t36 + t53 * t44;
t24 = qJDD(2) * pkin(2) + t59;
t54 = qJD(2) ^ 2;
t65 = t53 * t36 + t50 * t44;
t25 = -t54 * pkin(2) + t65;
t49 = sin(qJ(3));
t52 = cos(qJ(3));
t60 = t52 * t24 - t49 * t25;
t19 = t42 * pkin(3) + t60;
t66 = t49 * t24 + t52 * t25;
t20 = -t41 * pkin(3) + t66;
t45 = sin(pkin(9));
t47 = cos(pkin(9));
t57 = t47 * t19 - t45 * t20;
t70 = (t48 * t31 - t51 * t32) * t43 + m(6) * (-t42 * pkin(4) - t41 * pkin(7) - t57) - t28 * mrSges(6,1) + t27 * mrSges(6,2);
t67 = t45 * t19 + t47 * t20;
t10 = m(5) * t57 + t42 * mrSges(5,1) - t41 * mrSges(5,2) - t70;
t16 = -t41 * pkin(4) + t42 * pkin(7) + t67;
t26 = (-mrSges(6,1) * t51 + mrSges(6,2) * t48) * t43;
t35 = t46 * g(1) - t64 * g(2);
t33 = qJDD(4) - t35;
t13 = m(6) * (-t48 * t16 + t51 * t33) - t27 * mrSges(6,3) + qJDD(5) * mrSges(6,1) - t26 * t69 + qJD(5) * t32;
t14 = m(6) * (t51 * t16 + t48 * t33) + t28 * mrSges(6,3) - qJDD(5) * mrSges(6,2) + t26 * t68 - qJD(5) * t31;
t8 = m(5) * t67 - t41 * mrSges(5,1) - t42 * mrSges(5,2) - t48 * t13 + t51 * t14;
t6 = m(4) * t60 + t42 * mrSges(4,1) - t41 * mrSges(4,2) + t47 * t10 + t45 * t8;
t7 = m(4) * t66 - t41 * mrSges(4,1) - t42 * mrSges(4,2) - t45 * t10 + t47 * t8;
t4 = m(3) * t59 + qJDD(2) * mrSges(3,1) - t54 * mrSges(3,2) + t49 * t7 + t52 * t6;
t5 = m(3) * t65 - t54 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t49 * t6 + t52 * t7;
t62 = m(2) * t44 + t53 * t4 + t50 * t5;
t61 = m(5) * t33 + t51 * t13 + t48 * t14;
t58 = m(4) * t35 - t61;
t9 = (m(2) + m(3)) * t35 + t58;
t1 = m(2) * t36 - t50 * t4 + t53 * t5;
t2 = [-m(1) * g(1) + t64 * t1 - t46 * t9, t1, t5, t7, t8, t14; -m(1) * g(2) + t46 * t1 + t64 * t9, t9, t4, t6, t10, t13; -m(1) * g(3) + t62, t62, -m(3) * t35 - t58, -t58, t61, t70;];
f_new = t2;
