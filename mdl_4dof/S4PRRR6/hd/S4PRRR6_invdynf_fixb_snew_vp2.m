% Calculate vector of cutting forces with Newton-Euler
% S4PRRR6
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4PRRR6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR6_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR6_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR6_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR6_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR6_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:34:41
% EndTime: 2019-12-31 16:34:42
% DurationCPUTime: 0.37s
% Computational Cost: add. (2800->91), mult. (5509->128), div. (0->0), fcn. (3336->8), ass. (0->53)
t48 = sin(qJ(3));
t51 = cos(qJ(3));
t63 = qJD(2) * qJD(3);
t61 = t51 * t63;
t32 = t48 * qJDD(2) + t61;
t33 = t51 * qJDD(2) - t48 * t63;
t65 = qJD(2) * t48;
t36 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t65;
t64 = qJD(2) * t51;
t37 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t64;
t53 = qJD(2) ^ 2;
t47 = sin(qJ(4));
t50 = cos(qJ(4));
t26 = (t47 * t51 + t48 * t50) * qJD(2);
t16 = -t26 * qJD(4) - t47 * t32 + t50 * t33;
t25 = (-t47 * t48 + t50 * t51) * qJD(2);
t17 = t25 * qJD(4) + t50 * t32 + t47 * t33;
t43 = qJD(3) + qJD(4);
t23 = -t43 * mrSges(5,2) + t25 * mrSges(5,3);
t24 = t43 * mrSges(5,1) - t26 * mrSges(5,3);
t38 = qJD(3) * pkin(3) - pkin(6) * t65;
t44 = t51 ^ 2;
t46 = sin(pkin(7));
t66 = cos(pkin(7));
t35 = -t66 * g(1) - t46 * g(2);
t45 = -g(3) + qJDD(1);
t49 = sin(qJ(2));
t52 = cos(qJ(2));
t59 = -t49 * t35 + t52 * t45;
t56 = -qJDD(2) * pkin(2) - t59;
t55 = t16 * mrSges(5,1) + t25 * t23 - m(5) * (t38 * t65 - t33 * pkin(3) + (-pkin(6) * t44 - pkin(5)) * t53 + t56) - t17 * mrSges(5,2) - t26 * t24;
t69 = (t48 * t36 - t51 * t37) * qJD(2) + m(4) * (-t53 * pkin(5) + t56) - t33 * mrSges(4,1) + t32 * mrSges(4,2) - t55;
t67 = t52 * t35 + t49 * t45;
t22 = -t53 * pkin(2) + qJDD(2) * pkin(5) + t67;
t34 = t46 * g(1) - t66 * g(2);
t68 = t51 * t22 - t48 * t34;
t60 = -t48 * t22 - t51 * t34;
t11 = (-t32 + t61) * pkin(6) + (t48 * t51 * t53 + qJDD(3)) * pkin(3) + t60;
t12 = -t44 * t53 * pkin(3) + t33 * pkin(6) - qJD(3) * t38 + t68;
t19 = -t25 * mrSges(5,1) + t26 * mrSges(5,2);
t42 = qJDD(3) + qJDD(4);
t10 = m(5) * (t47 * t11 + t50 * t12) + t16 * mrSges(5,3) - t42 * mrSges(5,2) + t25 * t19 - t43 * t24;
t31 = (-mrSges(4,1) * t51 + mrSges(4,2) * t48) * qJD(2);
t9 = m(5) * (t50 * t11 - t47 * t12) - t17 * mrSges(5,3) + t42 * mrSges(5,1) - t26 * t19 + t43 * t23;
t5 = m(4) * t60 + qJDD(3) * mrSges(4,1) - t32 * mrSges(4,3) + qJD(3) * t37 + t47 * t10 - t31 * t65 + t50 * t9;
t6 = m(4) * t68 - qJDD(3) * mrSges(4,2) + t33 * mrSges(4,3) - qJD(3) * t36 + t50 * t10 + t31 * t64 - t47 * t9;
t3 = m(3) * t67 - t53 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t48 * t5 + t51 * t6;
t8 = m(3) * t59 + qJDD(2) * mrSges(3,1) - t53 * mrSges(3,2) - t69;
t62 = m(2) * t45 + t49 * t3 + t52 * t8;
t58 = -t48 * t6 - t51 * t5;
t4 = (m(2) + m(3)) * t34 + t58;
t1 = m(2) * t35 + t52 * t3 - t49 * t8;
t2 = [-m(1) * g(1) + t66 * t1 - t46 * t4, t1, t3, t6, t10; -m(1) * g(2) + t46 * t1 + t66 * t4, t4, t8, t5, t9; -m(1) * g(3) + t62, t62, -m(3) * t34 - t58, t69, -t55;];
f_new = t2;
