% Calculate vector of cutting forces with Newton-Euler
% S5PRRRR2
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
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
% Datum: 2019-12-05 17:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRRRR2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR2_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR2_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_invdynf_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR2_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR2_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR2_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:04:45
% EndTime: 2019-12-05 17:04:46
% DurationCPUTime: 0.36s
% Computational Cost: add. (4641->72), mult. (5229->98), div. (0->0), fcn. (2728->8), ass. (0->52)
t39 = qJDD(2) + qJDD(3);
t31 = qJDD(4) + t39;
t42 = sin(qJ(5));
t46 = cos(qJ(5));
t40 = qJD(2) + qJD(3);
t32 = qJD(4) + t40;
t62 = qJD(5) * t32;
t21 = t42 * t31 + t46 * t62;
t22 = t46 * t31 - t42 * t62;
t66 = t32 * t42;
t26 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t66;
t65 = t32 * t46;
t27 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t65;
t30 = t32 ^ 2;
t45 = sin(qJ(2));
t49 = cos(qJ(2));
t59 = t45 * g(1) - t49 * g(2);
t28 = qJDD(2) * pkin(2) + t59;
t50 = qJD(2) ^ 2;
t55 = -t49 * g(1) - t45 * g(2);
t29 = -t50 * pkin(2) + t55;
t44 = sin(qJ(3));
t48 = cos(qJ(3));
t58 = t48 * t28 - t44 * t29;
t18 = t39 * pkin(3) + t58;
t38 = t40 ^ 2;
t63 = t44 * t28 + t48 * t29;
t19 = -t38 * pkin(3) + t63;
t43 = sin(qJ(4));
t47 = cos(qJ(4));
t54 = t47 * t18 - t43 * t19;
t69 = (t42 * t26 - t46 * t27) * t32 + m(6) * (-t30 * pkin(6) - t54) - t22 * mrSges(6,1) + t21 * mrSges(6,2);
t68 = -m(1) - m(2);
t64 = t43 * t18 + t47 * t19;
t14 = t31 * pkin(6) + t64;
t20 = (-mrSges(6,1) * t46 + mrSges(6,2) * t42) * t32;
t41 = -g(3) + qJDD(1);
t12 = m(6) * (-t42 * t14 + t46 * t41) - t21 * mrSges(6,3) + qJDD(5) * mrSges(6,1) - t20 * t66 + qJD(5) * t27;
t13 = m(6) * (t46 * t14 + t42 * t41) + t22 * mrSges(6,3) - qJDD(5) * mrSges(6,2) + t20 * t65 - qJD(5) * t26;
t8 = m(5) * t64 - t30 * mrSges(5,1) - t31 * mrSges(5,2) - t42 * t12 + t46 * t13;
t9 = m(5) * t54 + t31 * mrSges(5,1) - t30 * mrSges(5,2) - t69;
t6 = m(4) * t58 + t39 * mrSges(4,1) - t38 * mrSges(4,2) + t43 * t8 + t47 * t9;
t7 = m(4) * t63 - t38 * mrSges(4,1) - t39 * mrSges(4,2) - t43 * t9 + t47 * t8;
t4 = m(3) * t59 + qJDD(2) * mrSges(3,1) - t50 * mrSges(3,2) + t44 * t7 + t48 * t6;
t5 = m(3) * t55 - t50 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t44 * t6 + t48 * t7;
t67 = t49 * t4 + t45 * t5;
t61 = m(5) * t41 + t46 * t12 + t42 * t13;
t60 = -t45 * t4 + t49 * t5;
t57 = m(4) * t41 + t61;
t56 = m(3) * t41 + t57;
t52 = m(2) * t41 + t56;
t1 = [t68 * g(1) + t60, -m(2) * g(1) + t60, t5, t7, t8, t13; t68 * g(2) + t67, -m(2) * g(2) + t67, t4, t6, t9, t12; -m(1) * g(3) + t52, t52, t56, t57, t61, t69;];
f_new = t1;
