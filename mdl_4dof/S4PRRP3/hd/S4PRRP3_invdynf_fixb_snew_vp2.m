% Calculate vector of cutting forces with Newton-Euler
% S4PRRP3
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
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-12-31 16:27
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4PRRP3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP3_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP3_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_invdynf_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP3_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP3_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP3_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:26:44
% EndTime: 2019-12-31 16:26:45
% DurationCPUTime: 0.28s
% Computational Cost: add. (1220->87), mult. (2397->113), div. (0->0), fcn. (1191->6), ass. (0->45)
t44 = sin(qJ(3));
t46 = cos(qJ(3));
t59 = qJD(2) * qJD(3);
t53 = t46 * t59;
t25 = t44 * qJDD(2) + t53;
t26 = t46 * qJDD(2) - t44 * t59;
t61 = qJD(2) * t44;
t32 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t61;
t60 = qJD(2) * t46;
t33 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t60;
t34 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t60;
t48 = qJD(2) ^ 2;
t42 = sin(pkin(6));
t43 = cos(pkin(6));
t28 = t42 * g(1) - t43 * g(2);
t29 = -t43 * g(1) - t42 * g(2);
t45 = sin(qJ(2));
t47 = cos(qJ(2));
t51 = t47 * t28 - t45 * t29;
t50 = -qJDD(2) * pkin(2) - t51;
t30 = qJD(3) * pkin(3) - qJ(4) * t61;
t31 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t61;
t40 = t46 ^ 2;
t54 = m(5) * (t30 * t61 - t26 * pkin(3) + qJDD(4) + (-qJ(4) * t40 - pkin(5)) * t48 + t50) + t31 * t61 + t25 * mrSges(5,2);
t70 = -(-t44 * t32 + (t33 + t34) * t46) * qJD(2) - (mrSges(4,1) + mrSges(5,1)) * t26 + m(4) * (-t48 * pkin(5) + t50) + t25 * mrSges(4,2) + t54;
t67 = pkin(3) * t48;
t63 = t45 * t28 + t47 * t29;
t15 = -t48 * pkin(2) + qJDD(2) * pkin(5) + t63;
t41 = -g(3) + qJDD(1);
t64 = t46 * t15 + t44 * t41;
t58 = qJD(2) * qJD(4);
t23 = (-mrSges(5,1) * t46 + mrSges(5,2) * t44) * qJD(2);
t24 = (-mrSges(4,1) * t46 + mrSges(4,2) * t44) * qJD(2);
t36 = t46 * t41;
t56 = qJD(3) * t33 + qJDD(3) * mrSges(5,1) + m(5) * (qJDD(3) * pkin(3) + t36 + (-t25 + t53) * qJ(4) + (t46 * t67 - t15 - 0.2e1 * t58) * t44);
t7 = m(4) * (-t44 * t15 + t36) + qJDD(3) * mrSges(4,1) + qJD(3) * t34 + (-mrSges(4,3) - mrSges(5,3)) * t25 + (-t23 - t24) * t61 + t56;
t55 = m(5) * (t26 * qJ(4) - qJD(3) * t30 - t40 * t67 + 0.2e1 * t46 * t58 + t64) + t23 * t60 + t26 * mrSges(5,3);
t8 = m(4) * t64 + t26 * mrSges(4,3) + t24 * t60 + (-mrSges(4,2) - mrSges(5,2)) * qJDD(3) + (-t32 - t31) * qJD(3) + t55;
t57 = m(3) * t41 + t44 * t8 + t46 * t7;
t52 = m(2) * t41 + t57;
t4 = m(3) * t51 + qJDD(2) * mrSges(3,1) - t48 * mrSges(3,2) - t70;
t3 = m(3) * t63 - t48 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t44 * t7 + t46 * t8;
t2 = m(2) * t29 + t47 * t3 - t45 * t4;
t1 = m(2) * t28 + t45 * t3 + t47 * t4;
t5 = [-m(1) * g(1) - t42 * t1 + t43 * t2, t2, t3, t8, -qJDD(3) * mrSges(5,2) - qJD(3) * t31 + t55; -m(1) * g(2) + t43 * t1 + t42 * t2, t1, t4, t7, -t25 * mrSges(5,3) - t23 * t61 + t56; -m(1) * g(3) + t52, t52, t57, t70, -t26 * mrSges(5,1) - t33 * t60 + t54;];
f_new = t5;
