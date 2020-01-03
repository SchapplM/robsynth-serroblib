% Calculate vector of cutting forces with Newton-Euler
% S4PRRR4
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
% Datum: 2019-12-31 16:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4PRRR4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR4_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR4_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR4_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR4_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR4_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR4_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR4_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:32:31
% EndTime: 2019-12-31 16:32:32
% DurationCPUTime: 0.36s
% Computational Cost: add. (2857->90), mult. (5677->128), div. (0->0), fcn. (3457->8), ass. (0->53)
t49 = sin(qJ(3));
t52 = cos(qJ(3));
t64 = qJD(2) * qJD(3);
t62 = t52 * t64;
t31 = t49 * qJDD(2) + t62;
t32 = t52 * qJDD(2) - t49 * t64;
t66 = qJD(2) * t49;
t35 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t66;
t65 = qJD(2) * t52;
t36 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t65;
t54 = qJD(2) ^ 2;
t48 = sin(qJ(4));
t51 = cos(qJ(4));
t26 = (t48 * t52 + t49 * t51) * qJD(2);
t16 = -t26 * qJD(4) - t48 * t31 + t51 * t32;
t25 = (-t48 * t49 + t51 * t52) * qJD(2);
t17 = t25 * qJD(4) + t51 * t31 + t48 * t32;
t43 = qJD(3) + qJD(4);
t23 = -t43 * mrSges(5,2) + t25 * mrSges(5,3);
t24 = t43 * mrSges(5,1) - t26 * mrSges(5,3);
t37 = qJD(3) * pkin(3) - pkin(6) * t66;
t44 = t52 ^ 2;
t46 = sin(pkin(7));
t47 = cos(pkin(7));
t33 = t46 * g(1) - t47 * g(2);
t34 = -t47 * g(1) - t46 * g(2);
t50 = sin(qJ(2));
t53 = cos(qJ(2));
t59 = t53 * t33 - t50 * t34;
t57 = -qJDD(2) * pkin(2) - t59;
t56 = t16 * mrSges(5,1) + t25 * t23 - m(5) * (t37 * t66 - t32 * pkin(3) + (-pkin(6) * t44 - pkin(5)) * t54 + t57) - t17 * mrSges(5,2) - t26 * t24;
t69 = (t49 * t35 - t52 * t36) * qJD(2) + m(4) * (-t54 * pkin(5) + t57) - t32 * mrSges(4,1) + t31 * mrSges(4,2) - t56;
t67 = t50 * t33 + t53 * t34;
t22 = -t54 * pkin(2) + qJDD(2) * pkin(5) + t67;
t45 = -g(3) + qJDD(1);
t68 = t52 * t22 + t49 * t45;
t60 = -t49 * t22 + t52 * t45;
t11 = (-t31 + t62) * pkin(6) + (t49 * t52 * t54 + qJDD(3)) * pkin(3) + t60;
t12 = -t44 * t54 * pkin(3) + t32 * pkin(6) - qJD(3) * t37 + t68;
t20 = -t25 * mrSges(5,1) + t26 * mrSges(5,2);
t42 = qJDD(3) + qJDD(4);
t10 = m(5) * (t48 * t11 + t51 * t12) + t16 * mrSges(5,3) - t42 * mrSges(5,2) + t25 * t20 - t43 * t24;
t30 = (-mrSges(4,1) * t52 + mrSges(4,2) * t49) * qJD(2);
t9 = m(5) * (t51 * t11 - t48 * t12) - t17 * mrSges(5,3) + t42 * mrSges(5,1) - t26 * t20 + t43 * t23;
t6 = m(4) * t60 + qJDD(3) * mrSges(4,1) - t31 * mrSges(4,3) + qJD(3) * t36 + t48 * t10 - t30 * t66 + t51 * t9;
t7 = m(4) * t68 - qJDD(3) * mrSges(4,2) + t32 * mrSges(4,3) - qJD(3) * t35 + t51 * t10 + t30 * t65 - t48 * t9;
t63 = m(3) * t45 + t49 * t7 + t52 * t6;
t61 = m(2) * t45 + t63;
t8 = m(3) * t59 + qJDD(2) * mrSges(3,1) - t54 * mrSges(3,2) - t69;
t3 = m(3) * t67 - t54 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t49 * t6 + t52 * t7;
t2 = m(2) * t34 + t53 * t3 - t50 * t8;
t1 = m(2) * t33 + t50 * t3 + t53 * t8;
t4 = [-m(1) * g(1) - t46 * t1 + t47 * t2, t2, t3, t7, t10; -m(1) * g(2) + t47 * t1 + t46 * t2, t1, t8, t6, t9; -m(1) * g(3) + t61, t61, t63, t69, -t56;];
f_new = t4;
