% Calculate vector of cutting forces with Newton-Euler
% S4RPRR8
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% f_new [3x5]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4RPRR8_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR8_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR8_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR8_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR8_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR8_invdynf_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR8_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR8_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR8_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:55:06
% EndTime: 2019-12-31 16:55:07
% DurationCPUTime: 0.32s
% Computational Cost: add. (1993->97), mult. (3856->127), div. (0->0), fcn. (2076->6), ass. (0->52)
t46 = sin(qJ(1));
t49 = cos(qJ(1));
t58 = -t49 * g(1) - t46 * g(2);
t56 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t58;
t71 = -m(2) - m(3);
t70 = -pkin(1) - pkin(5);
t69 = (mrSges(2,1) - mrSges(3,2));
t68 = -mrSges(2,2) + mrSges(3,3);
t50 = qJD(1) ^ 2;
t61 = t46 * g(1) - t49 * g(2);
t54 = -t50 * qJ(2) + qJDD(2) - t61;
t23 = t70 * qJDD(1) + t54;
t45 = sin(qJ(3));
t48 = cos(qJ(3));
t67 = t45 * g(3) + t48 * t23;
t66 = qJD(1) * t45;
t65 = qJD(1) * t48;
t64 = qJD(1) * qJD(3);
t31 = (mrSges(4,1) * t45 + mrSges(4,2) * t48) * qJD(1);
t59 = t45 * t64;
t33 = t48 * qJDD(1) - t59;
t34 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t66;
t44 = sin(qJ(4));
t47 = cos(qJ(4));
t26 = (-t44 * t48 - t45 * t47) * qJD(1);
t32 = -t45 * qJDD(1) - t48 * t64;
t14 = t26 * qJD(4) + t44 * t32 + t47 * t33;
t27 = (-t44 * t45 + t47 * t48) * qJD(1);
t16 = -t26 * mrSges(5,1) + t27 * mrSges(5,2);
t41 = qJD(3) + qJD(4);
t20 = -t41 * mrSges(5,2) + t26 * mrSges(5,3);
t40 = qJDD(3) + qJDD(4);
t8 = (-t33 - t59) * pkin(6) + (-t45 * t48 * t50 + qJDD(3)) * pkin(3) + t67;
t36 = qJD(3) * pkin(3) - pkin(6) * t65;
t43 = t45 ^ 2;
t60 = -t48 * g(3) + t45 * t23;
t9 = -t43 * t50 * pkin(3) + t32 * pkin(6) - qJD(3) * t36 + t60;
t6 = m(5) * (-t44 * t9 + t47 * t8) - t14 * mrSges(5,3) + t40 * mrSges(5,1) - t27 * t16 + t41 * t20;
t13 = -t27 * qJD(4) + t47 * t32 - t44 * t33;
t21 = t41 * mrSges(5,1) - t27 * mrSges(5,3);
t7 = m(5) * (t44 * t8 + t47 * t9) + t13 * mrSges(5,3) - t40 * mrSges(5,2) + t26 * t16 - t41 * t21;
t3 = m(4) * t67 + qJDD(3) * mrSges(4,1) - t33 * mrSges(4,3) + qJD(3) * t34 - t31 * t65 + t44 * t7 + t47 * t6;
t35 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t65;
t4 = m(4) * t60 - qJDD(3) * mrSges(4,2) + t32 * mrSges(4,3) - qJD(3) * t35 - t31 * t66 - t44 * t6 + t47 * t7;
t62 = -t45 * t3 + t48 * t4;
t55 = -m(3) * (-qJDD(1) * pkin(1) + t54) - t48 * t3 - t45 * t4;
t53 = -t13 * mrSges(5,1) - t26 * t20 + m(5) * (t36 * t65 - t32 * pkin(3) + (-pkin(6) * t43 + t70) * t50 + t56) + t14 * mrSges(5,2) + t27 * t21;
t52 = -t32 * mrSges(4,1) + m(4) * (t70 * t50 + t56) + t34 * t66 + t35 * t65 + t33 * mrSges(4,2) + t53;
t51 = -m(3) * (t50 * pkin(1) - t56) + t52;
t5 = m(2) * t58 + t68 * qJDD(1) - (t69 * t50) + t51;
t1 = m(2) * t61 + t69 * qJDD(1) + t68 * t50 + t55;
t2 = [-m(1) * g(1) - t46 * t1 + t49 * t5, t5, -m(3) * g(3) + t62, t4, t7; -m(1) * g(2) + t49 * t1 + t46 * t5, t1, -(t50 * mrSges(3,2)) - qJDD(1) * mrSges(3,3) - t51, t3, t6; (-m(1) + t71) * g(3) + t62, t71 * g(3) + t62, qJDD(1) * mrSges(3,2) - t50 * mrSges(3,3) - t55, t52, t53;];
f_new = t2;
