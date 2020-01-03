% Calculate vector of cutting forces with Newton-Euler
% S4RPRP7
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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
% Datum: 2019-12-31 16:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4RPRP7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP7_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP7_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP7_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP7_invdynf_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP7_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP7_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP7_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:47:03
% EndTime: 2019-12-31 16:47:04
% DurationCPUTime: 0.24s
% Computational Cost: add. (835->92), mult. (1548->115), div. (0->0), fcn. (592->4), ass. (0->43)
t36 = sin(qJ(1));
t38 = cos(qJ(1));
t48 = -t38 * g(1) - t36 * g(2);
t65 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t48;
t40 = qJD(1) ^ 2;
t51 = t36 * g(1) - t38 * g(2);
t45 = -t40 * qJ(2) + qJDD(2) - t51;
t62 = -pkin(1) - pkin(5);
t13 = t62 * qJDD(1) + t45;
t35 = sin(qJ(3));
t37 = cos(qJ(3));
t21 = (pkin(3) * t35 - qJ(4) * t37) * qJD(1);
t39 = qJD(3) ^ 2;
t61 = t35 * g(3);
t64 = m(5) * (-qJDD(3) * pkin(3) - t61 - t39 * qJ(4) + qJDD(4) + (qJD(1) * t21 - t13) * t37);
t63 = -m(2) - m(3);
t60 = (mrSges(2,1) - mrSges(3,2));
t59 = -mrSges(2,2) + mrSges(3,3);
t58 = -mrSges(4,3) - mrSges(5,2);
t57 = qJD(1) * t35;
t56 = qJD(1) * t37;
t55 = qJD(1) * qJD(3);
t29 = -(qJD(3) * mrSges(5,1)) + mrSges(5,2) * t56;
t50 = -t37 * g(3) + t35 * t13;
t53 = qJD(3) * t29 + qJDD(3) * mrSges(5,3) + m(5) * (-t39 * pkin(3) + qJDD(3) * qJ(4) + (2 * qJD(4) * qJD(3)) - t21 * t57 + t50);
t24 = t35 * qJDD(1) + t37 * t55;
t28 = (qJD(3) * mrSges(4,1)) - mrSges(4,3) * t56;
t22 = (mrSges(5,1) * t35 - mrSges(5,3) * t37) * qJD(1);
t49 = qJD(1) * (-t22 - (mrSges(4,1) * t35 + mrSges(4,2) * t37) * qJD(1));
t4 = m(4) * t50 - qJDD(3) * mrSges(4,2) - qJD(3) * t28 + t58 * t24 + t35 * t49 + t53;
t25 = t37 * qJDD(1) - t35 * t55;
t27 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t57;
t30 = -mrSges(5,2) * t57 + qJD(3) * mrSges(5,3);
t5 = m(4) * (t37 * t13 + t61) - t64 + t58 * t25 + (mrSges(4,1) + mrSges(5,1)) * qJDD(3) + (t27 + t30) * qJD(3) + t37 * t49;
t52 = -t35 * t5 + t37 * t4;
t46 = -m(3) * (-qJDD(1) * pkin(1) + t45) - t35 * t4 - t37 * t5;
t43 = t62 * t40 - t65;
t44 = -t25 * mrSges(5,3) - t29 * t56 + t30 * t57 + t24 * mrSges(5,1) + m(5) * (t24 * pkin(3) - t25 * qJ(4) + (-0.2e1 * qJD(4) * t37 + (pkin(3) * t37 + qJ(4) * t35) * qJD(3)) * qJD(1) + t43);
t42 = m(4) * t43 + t24 * mrSges(4,1) + t25 * mrSges(4,2) + t27 * t57 + t28 * t56 + t44;
t41 = -m(3) * (t40 * pkin(1) + t65) + t42;
t2 = m(2) * t48 + t59 * qJDD(1) - (t60 * t40) + t41;
t1 = m(2) * t51 + t60 * qJDD(1) + t59 * t40 + t46;
t3 = [-m(1) * g(1) - t36 * t1 + t38 * t2, t2, -m(3) * g(3) + t52, t4, -t24 * mrSges(5,2) - t22 * t57 + t53; -m(1) * g(2) + t38 * t1 + t36 * t2, t1, -(t40 * mrSges(3,2)) - qJDD(1) * mrSges(3,3) - t41, t5, t44; (-m(1) + t63) * g(3) + t52, t63 * g(3) + t52, qJDD(1) * mrSges(3,2) - t40 * mrSges(3,3) - t46, t42, -qJDD(3) * mrSges(5,1) + t25 * mrSges(5,2) - qJD(3) * t30 + t22 * t56 + t64;];
f_new = t3;
