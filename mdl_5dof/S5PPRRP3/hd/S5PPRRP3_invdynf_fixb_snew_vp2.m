% Calculate vector of cutting forces with Newton-Euler
% S5PPRRP3
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PPRRP3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP3_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP3_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP3_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP3_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP3_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:10:36
% EndTime: 2019-12-05 15:10:36
% DurationCPUTime: 0.39s
% Computational Cost: add. (2563->95), mult. (4409->127), div. (0->0), fcn. (2513->8), ass. (0->51)
t51 = qJD(3) ^ 2;
t43 = sin(pkin(7));
t45 = cos(pkin(7));
t33 = -t45 * g(1) - t43 * g(2);
t41 = -g(3) + qJDD(1);
t42 = sin(pkin(8));
t44 = cos(pkin(8));
t22 = t44 * t33 + t42 * t41;
t55 = t43 * g(1) - t45 * g(2);
t32 = qJDD(2) - t55;
t47 = sin(qJ(3));
t49 = cos(qJ(3));
t57 = -t47 * t22 + t49 * t32;
t16 = -qJDD(3) * pkin(3) - t51 * pkin(6) - t57;
t46 = sin(qJ(4));
t48 = cos(qJ(4));
t59 = qJD(3) * qJD(4);
t29 = t46 * qJDD(3) + t48 * t59;
t30 = t48 * qJDD(3) - t46 * t59;
t61 = qJD(3) * t46;
t34 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t61;
t35 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t61;
t60 = qJD(3) * t48;
t37 = mrSges(6,2) * t60 + qJD(4) * mrSges(6,3);
t62 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t60 + t37;
t66 = m(6) * (-t30 * pkin(4) - t29 * qJ(5) + (-0.2e1 * qJD(5) * t46 + (pkin(4) * t46 - qJ(5) * t48) * qJD(4)) * qJD(3) + t16) - t30 * mrSges(6,1);
t73 = (-t62 * t48 + (t34 - t35) * t46) * qJD(3) + m(5) * t16 - t30 * mrSges(5,1) + (mrSges(5,2) - mrSges(6,3)) * t29 + t66;
t64 = t49 * t22 + t47 * t32;
t17 = -t51 * pkin(3) + qJDD(3) * pkin(6) + t64;
t26 = (-pkin(4) * t48 - qJ(5) * t46) * qJD(3);
t50 = qJD(4) ^ 2;
t21 = t42 * t33 - t44 * t41;
t69 = t48 * t21;
t70 = m(6) * (-qJDD(4) * pkin(4) - t50 * qJ(5) - t69 + qJDD(5) + (qJD(3) * t26 + t17) * t46);
t67 = mrSges(5,3) + mrSges(6,2);
t65 = t48 * t17 + t46 * t21;
t27 = (-mrSges(6,1) * t48 - mrSges(6,3) * t46) * qJD(3);
t28 = (-mrSges(5,1) * t48 + mrSges(5,2) * t46) * qJD(3);
t10 = m(5) * (-t46 * t17 + t69) - t70 - t67 * t29 + (mrSges(5,1) + mrSges(6,1)) * qJDD(4) + t62 * qJD(4) + (-t27 - t28) * t61;
t56 = m(6) * (-t50 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t26 * t60 + t65) + t27 * t60 + qJD(4) * t35 + qJDD(4) * mrSges(6,3);
t9 = m(5) * t65 - qJDD(4) * mrSges(5,2) - qJD(4) * t34 + t28 * t60 + t67 * t30 + t56;
t6 = m(4) * t64 - t51 * mrSges(4,1) - qJDD(3) * mrSges(4,2) - t46 * t10 + t48 * t9;
t8 = m(4) * t57 + qJDD(3) * mrSges(4,1) - t51 * mrSges(4,2) - t73;
t4 = m(3) * t22 - t47 * t8 + t49 * t6;
t54 = -t48 * t10 - t46 * t9;
t7 = (-m(3) - m(4)) * t21 + t54;
t58 = m(2) * t41 + t42 * t4 + t44 * t7;
t52 = m(3) * t32 + t47 * t6 + t49 * t8;
t3 = m(2) * t55 - t52;
t1 = m(2) * t33 + t44 * t4 - t42 * t7;
t2 = [-m(1) * g(1) + t45 * t1 - t43 * t3, t1, t4, t6, t9, t30 * mrSges(6,2) + t56; -m(1) * g(2) + t43 * t1 + t45 * t3, t3, t7, t8, t10, -t29 * mrSges(6,3) + (-t46 * t35 - t48 * t37) * qJD(3) + t66; -m(1) * g(3) + t58, t58, t52, m(4) * t21 - t54, t73, -qJDD(4) * mrSges(6,1) + t29 * mrSges(6,2) - qJD(4) * t37 + t27 * t61 + t70;];
f_new = t2;
