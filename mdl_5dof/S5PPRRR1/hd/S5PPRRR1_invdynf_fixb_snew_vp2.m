% Calculate vector of cutting forces with Newton-Euler
% S5PPRRR1
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
% Datum: 2019-12-05 15:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PPRRR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR1_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR1_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR1_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR1_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR1_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:12:41
% EndTime: 2019-12-05 15:12:41
% DurationCPUTime: 0.45s
% Computational Cost: add. (5190->70), mult. (6794->101), div. (0->0), fcn. (4710->10), ass. (0->51)
t37 = qJDD(3) + qJDD(4);
t44 = sin(qJ(5));
t47 = cos(qJ(5));
t38 = qJD(3) + qJD(4);
t60 = qJD(5) * t38;
t27 = t44 * t37 + t47 * t60;
t28 = t47 * t37 - t44 * t60;
t64 = t38 * t44;
t31 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t64;
t63 = t38 * t47;
t32 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t63;
t36 = t38 ^ 2;
t41 = sin(pkin(8));
t43 = cos(pkin(8));
t34 = -t43 * g(1) - t41 * g(2);
t39 = -g(3) + qJDD(1);
t40 = sin(pkin(9));
t42 = cos(pkin(9));
t24 = -t40 * t34 + t42 * t39;
t25 = t42 * t34 + t40 * t39;
t46 = sin(qJ(3));
t49 = cos(qJ(3));
t57 = t49 * t24 - t46 * t25;
t19 = qJDD(3) * pkin(3) + t57;
t50 = qJD(3) ^ 2;
t61 = t46 * t24 + t49 * t25;
t20 = -t50 * pkin(3) + t61;
t45 = sin(qJ(4));
t48 = cos(qJ(4));
t54 = t48 * t19 - t45 * t20;
t65 = (t44 * t31 - t47 * t32) * t38 + m(6) * (-t37 * pkin(4) - t36 * pkin(7) - t54) - t28 * mrSges(6,1) + t27 * mrSges(6,2);
t62 = t45 * t19 + t48 * t20;
t10 = m(5) * t54 + t37 * mrSges(5,1) - t36 * mrSges(5,2) - t65;
t16 = -t36 * pkin(4) + t37 * pkin(7) + t62;
t26 = (-mrSges(6,1) * t47 + mrSges(6,2) * t44) * t38;
t55 = t41 * g(1) - t43 * g(2);
t33 = qJDD(2) - t55;
t13 = m(6) * (-t44 * t16 + t47 * t33) - t27 * mrSges(6,3) + qJDD(5) * mrSges(6,1) - t26 * t64 + qJD(5) * t32;
t14 = m(6) * (t47 * t16 + t44 * t33) + t28 * mrSges(6,3) - qJDD(5) * mrSges(6,2) + t26 * t63 - qJD(5) * t31;
t8 = m(5) * t62 - t36 * mrSges(5,1) - t37 * mrSges(5,2) - t44 * t13 + t47 * t14;
t6 = m(4) * t57 + qJDD(3) * mrSges(4,1) - t50 * mrSges(4,2) + t48 * t10 + t45 * t8;
t7 = m(4) * t61 - t50 * mrSges(4,1) - qJDD(3) * mrSges(4,2) - t45 * t10 + t48 * t8;
t4 = m(3) * t24 + t46 * t7 + t49 * t6;
t5 = m(3) * t25 - t46 * t6 + t49 * t7;
t59 = m(2) * t39 + t42 * t4 + t40 * t5;
t58 = m(5) * t33 + t47 * t13 + t44 * t14;
t56 = m(4) * t33 + t58;
t52 = m(3) * t33 + t56;
t9 = m(2) * t55 - t52;
t1 = m(2) * t34 - t40 * t4 + t42 * t5;
t2 = [-m(1) * g(1) + t43 * t1 - t41 * t9, t1, t5, t7, t8, t14; -m(1) * g(2) + t41 * t1 + t43 * t9, t9, t4, t6, t10, t13; -m(1) * g(3) + t59, t59, t52, t56, t58, t65;];
f_new = t2;
