% Calculate vector of cutting forces with Newton-Euler
% S4PRRP5
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
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4PRRP5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP5_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP5_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_invdynf_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP5_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP5_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP5_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:28:52
% EndTime: 2019-12-31 16:28:52
% DurationCPUTime: 0.28s
% Computational Cost: add. (1176->88), mult. (2261->113), div. (0->0), fcn. (1098->6), ass. (0->45)
t43 = sin(qJ(3));
t45 = cos(qJ(3));
t58 = qJD(2) * qJD(3);
t52 = t45 * t58;
t26 = t43 * qJDD(2) + t52;
t27 = t45 * qJDD(2) - t43 * t58;
t60 = qJD(2) * t43;
t33 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t60;
t59 = qJD(2) * t45;
t34 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t59;
t35 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t59;
t47 = qJD(2) ^ 2;
t42 = sin(pkin(6));
t61 = cos(pkin(6));
t30 = -t61 * g(1) - t42 * g(2);
t41 = -g(3) + qJDD(1);
t44 = sin(qJ(2));
t46 = cos(qJ(2));
t51 = -t44 * t30 + t46 * t41;
t49 = -qJDD(2) * pkin(2) - t51;
t31 = qJD(3) * pkin(3) - qJ(4) * t60;
t32 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t60;
t40 = t45 ^ 2;
t53 = m(5) * (t31 * t60 - t27 * pkin(3) + qJDD(4) + (-qJ(4) * t40 - pkin(5)) * t47 + t49) + t32 * t60 + t26 * mrSges(5,2);
t70 = -(-t43 * t33 + (t34 + t35) * t45) * qJD(2) - (mrSges(4,1) + mrSges(5,1)) * t27 + m(4) * (-t47 * pkin(5) + t49) + t26 * mrSges(4,2) + t53;
t67 = pkin(3) * t47;
t63 = t46 * t30 + t44 * t41;
t15 = -t47 * pkin(2) + qJDD(2) * pkin(5) + t63;
t29 = t42 * g(1) - t61 * g(2);
t64 = t45 * t15 - t43 * t29;
t57 = qJD(2) * qJD(4);
t22 = t45 * t29;
t24 = (-mrSges(5,1) * t45 + mrSges(5,2) * t43) * qJD(2);
t25 = (-mrSges(4,1) * t45 + mrSges(4,2) * t43) * qJD(2);
t55 = qJD(3) * t34 + qJDD(3) * mrSges(5,1) + m(5) * (qJDD(3) * pkin(3) - t22 + (-t26 + t52) * qJ(4) + (t45 * t67 - t15 - 0.2e1 * t57) * t43);
t7 = m(4) * (-t43 * t15 - t22) + qJDD(3) * mrSges(4,1) + qJD(3) * t35 + (-mrSges(4,3) - mrSges(5,3)) * t26 + (-t24 - t25) * t60 + t55;
t54 = m(5) * (t27 * qJ(4) - qJD(3) * t31 - t40 * t67 + 0.2e1 * t45 * t57 + t64) + t24 * t59 + t27 * mrSges(5,3);
t8 = m(4) * t64 + t27 * mrSges(4,3) + t25 * t59 + (-mrSges(4,2) - mrSges(5,2)) * qJDD(3) + (-t33 - t32) * qJD(3) + t54;
t3 = m(3) * t63 - t47 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t43 * t7 + t45 * t8;
t6 = m(3) * t51 + qJDD(2) * mrSges(3,1) - t47 * mrSges(3,2) - t70;
t56 = m(2) * t41 + t44 * t3 + t46 * t6;
t50 = -t43 * t8 - t45 * t7;
t4 = (m(2) + m(3)) * t29 + t50;
t1 = m(2) * t30 + t46 * t3 - t44 * t6;
t2 = [-m(1) * g(1) + t61 * t1 - t42 * t4, t1, t3, t8, -qJDD(3) * mrSges(5,2) - qJD(3) * t32 + t54; -m(1) * g(2) + t42 * t1 + t61 * t4, t4, t6, t7, -t26 * mrSges(5,3) - t24 * t60 + t55; -m(1) * g(3) + t56, t56, -m(3) * t29 - t50, t70, -t27 * mrSges(5,1) - t34 * t59 + t53;];
f_new = t2;
