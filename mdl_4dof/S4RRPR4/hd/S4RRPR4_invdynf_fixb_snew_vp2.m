% Calculate vector of cutting forces with Newton-Euler
% S4RRPR4
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4RRPR4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR4_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR4_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR4_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR4_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR4_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:02:27
% EndTime: 2019-12-31 17:02:28
% DurationCPUTime: 0.36s
% Computational Cost: add. (3904->87), mult. (5435->116), div. (0->0), fcn. (3165->8), ass. (0->54)
t40 = qJD(1) + qJD(2);
t36 = t40 ^ 2;
t42 = cos(pkin(7));
t39 = t42 ^ 2;
t41 = sin(pkin(7));
t63 = t41 ^ 2 + t39;
t69 = t63 * mrSges(4,3);
t68 = -m(2) - m(3);
t37 = qJDD(1) + qJDD(2);
t45 = sin(qJ(1));
t48 = cos(qJ(1));
t61 = t45 * g(1) - t48 * g(2);
t31 = qJDD(1) * pkin(1) + t61;
t49 = qJD(1) ^ 2;
t57 = -t48 * g(1) - t45 * g(2);
t32 = -t49 * pkin(1) + t57;
t44 = sin(qJ(2));
t47 = cos(qJ(2));
t64 = t44 * t31 + t47 * t32;
t22 = -t36 * pkin(2) + t37 * qJ(3) + t64;
t62 = qJD(3) * t40;
t60 = -t42 * g(3) - 0.2e1 * t41 * t62;
t65 = pkin(6) * t37;
t66 = pkin(3) * t42;
t11 = (t36 * t66 - t22 - t65) * t41 + t60;
t58 = -t41 * g(3) + (t22 + 0.2e1 * t62) * t42;
t12 = -t39 * t36 * pkin(3) + t42 * t65 + t58;
t43 = sin(qJ(4));
t46 = cos(qJ(4));
t52 = -t41 * t43 + t42 * t46;
t25 = t52 * t40;
t53 = t41 * t46 + t42 * t43;
t26 = t53 * t40;
t18 = -t25 * mrSges(5,1) + t26 * mrSges(5,2);
t20 = -t26 * qJD(4) + t52 * t37;
t24 = qJD(4) * mrSges(5,1) - t26 * mrSges(5,3);
t10 = m(5) * (t43 * t11 + t46 * t12) + t20 * mrSges(5,3) - qJDD(4) * mrSges(5,2) + t25 * t18 - qJD(4) * t24;
t55 = -t42 * mrSges(4,1) + t41 * mrSges(4,2);
t54 = t37 * mrSges(4,3) + t36 * t55;
t21 = t25 * qJD(4) + t53 * t37;
t23 = -qJD(4) * mrSges(5,2) + t25 * mrSges(5,3);
t9 = m(5) * (t46 * t11 - t43 * t12) - t21 * mrSges(5,3) + qJDD(4) * mrSges(5,1) - t26 * t18 + qJD(4) * t23;
t6 = m(4) * t60 + t43 * t10 + t46 * t9 + (-m(4) * t22 - t54) * t41;
t7 = m(4) * t58 + t46 * t10 + t54 * t42 - t43 * t9;
t67 = t41 * t7 + t42 * t6;
t59 = t47 * t31 - t44 * t32;
t56 = qJDD(3) - t59;
t51 = t20 * mrSges(5,1) + t25 * t23 - m(5) * ((-pkin(2) - t66) * t37 + (-t63 * pkin(6) - qJ(3)) * t36 + t56) - t26 * t24 - t21 * mrSges(5,2);
t50 = m(4) * (-t37 * pkin(2) - t36 * qJ(3) + t56) - t51;
t8 = m(3) * t59 + (mrSges(3,1) - t55) * t37 + (-mrSges(3,2) + t69) * t36 - t50;
t3 = m(3) * t64 - t36 * mrSges(3,1) - t37 * mrSges(3,2) - t41 * t6 + t42 * t7;
t2 = m(2) * t57 - t49 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t47 * t3 - t44 * t8;
t1 = m(2) * t61 + qJDD(1) * mrSges(2,1) - t49 * mrSges(2,2) + t44 * t3 + t47 * t8;
t4 = [-m(1) * g(1) - t45 * t1 + t48 * t2, t2, t3, t7, t10; -m(1) * g(2) + t48 * t1 + t45 * t2, t1, t8, t6, t9; (-m(1) + t68) * g(3) + t67, t68 * g(3) + t67, -m(3) * g(3) + t67, -t36 * t69 + t55 * t37 + t50, -t51;];
f_new = t4;
