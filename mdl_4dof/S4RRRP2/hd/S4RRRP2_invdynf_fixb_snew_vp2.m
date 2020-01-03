% Calculate vector of cutting forces with Newton-Euler
% S4RRRP2
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
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4RRRP2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP2_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP2_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_invdynf_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP2_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP2_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP2_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:12:54
% EndTime: 2019-12-31 17:12:55
% DurationCPUTime: 0.32s
% Computational Cost: add. (2066->96), mult. (2646->120), div. (0->0), fcn. (1191->6), ass. (0->49)
t37 = qJDD(1) + qJDD(2);
t41 = sin(qJ(3));
t44 = cos(qJ(3));
t38 = qJD(1) + qJD(2);
t59 = qJD(3) * t38;
t54 = t44 * t59;
t22 = t41 * t37 + t54;
t23 = t44 * t37 - t41 * t59;
t65 = t38 * t41;
t32 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t65;
t64 = t38 * t44;
t33 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t64;
t34 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t64;
t36 = t38 ^ 2;
t43 = sin(qJ(1));
t46 = cos(qJ(1));
t53 = t43 * g(1) - t46 * g(2);
t28 = qJDD(1) * pkin(1) + t53;
t47 = qJD(1) ^ 2;
t50 = -t46 * g(1) - t43 * g(2);
t29 = -t47 * pkin(1) + t50;
t42 = sin(qJ(2));
t45 = cos(qJ(2));
t51 = t45 * t28 - t42 * t29;
t49 = -t37 * pkin(2) - t51;
t30 = qJD(3) * pkin(3) - qJ(4) * t65;
t31 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t65;
t40 = t44 ^ 2;
t55 = m(5) * (t30 * t65 - t23 * pkin(3) + qJDD(4) + (-qJ(4) * t40 - pkin(6)) * t36 + t49) + t31 * t65 + t22 * mrSges(5,2);
t72 = -(-t41 * t32 + (t33 + t34) * t44) * t38 - (mrSges(4,1) + mrSges(5,1)) * t23 + m(4) * (-t36 * pkin(6) + t49) + t22 * mrSges(4,2) + t55;
t69 = -m(2) - m(3);
t61 = t42 * t28 + t45 * t29;
t15 = -t36 * pkin(2) + t37 * pkin(6) + t61;
t20 = (-mrSges(5,1) * t44 + mrSges(5,2) * t41) * t38;
t21 = (-mrSges(4,1) * t44 + mrSges(4,2) * t41) * t38;
t58 = qJD(4) * t38;
t66 = t44 * g(3);
t67 = pkin(3) * t36;
t57 = qJD(3) * t33 + qJDD(3) * mrSges(5,1) + m(5) * (qJDD(3) * pkin(3) - t66 + (-t22 + t54) * qJ(4) + (t44 * t67 - t15 - 0.2e1 * t58) * t41);
t7 = m(4) * (-t41 * t15 - t66) + qJDD(3) * mrSges(4,1) + qJD(3) * t34 + (-t20 - t21) * t65 + (-mrSges(4,3) - mrSges(5,3)) * t22 + t57;
t52 = -t41 * g(3) + t44 * t15;
t56 = m(5) * (t23 * qJ(4) - qJD(3) * t30 - t40 * t67 + 0.2e1 * t44 * t58 + t52) + t20 * t64 + t23 * mrSges(5,3);
t8 = m(4) * t52 + t23 * mrSges(4,3) + t21 * t64 + (-mrSges(4,2) - mrSges(5,2)) * qJDD(3) + (-t32 - t31) * qJD(3) + t56;
t68 = t41 * t8 + t44 * t7;
t4 = m(3) * t51 + t37 * mrSges(3,1) - t36 * mrSges(3,2) - t72;
t3 = m(3) * t61 - t36 * mrSges(3,1) - t37 * mrSges(3,2) - t41 * t7 + t44 * t8;
t2 = m(2) * t50 - t47 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t45 * t3 - t42 * t4;
t1 = m(2) * t53 + qJDD(1) * mrSges(2,1) - t47 * mrSges(2,2) + t42 * t3 + t45 * t4;
t5 = [-m(1) * g(1) - t43 * t1 + t46 * t2, t2, t3, t8, -qJDD(3) * mrSges(5,2) - qJD(3) * t31 + t56; -m(1) * g(2) + t46 * t1 + t43 * t2, t1, t4, t7, -t22 * mrSges(5,3) - t20 * t65 + t57; (-m(1) + t69) * g(3) + t68, g(3) * t69 + t68, -m(3) * g(3) + t68, t72, -t23 * mrSges(5,1) - t33 * t64 + t55;];
f_new = t5;
