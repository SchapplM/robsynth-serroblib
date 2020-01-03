% Calculate vector of cutting forces with Newton-Euler
% S4PRPR3
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
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-31 16:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4PRPR3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR3_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR3_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR3_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR3_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR3_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR3_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:20:49
% EndTime: 2019-12-31 16:20:50
% DurationCPUTime: 0.32s
% Computational Cost: add. (2290->78), mult. (4982->109), div. (0->0), fcn. (3165->8), ass. (0->50)
t50 = qJD(2) ^ 2;
t44 = cos(pkin(7));
t40 = t44 ^ 2;
t42 = sin(pkin(7));
t64 = t42 ^ 2 + t40;
t68 = t64 * mrSges(4,3);
t67 = pkin(3) * t50;
t43 = sin(pkin(6));
t45 = cos(pkin(6));
t31 = t43 * g(1) - t45 * g(2);
t32 = -t45 * g(1) - t43 * g(2);
t47 = sin(qJ(2));
t49 = cos(qJ(2));
t66 = t47 * t31 + t49 * t32;
t41 = -g(3) + qJDD(1);
t62 = qJD(2) * qJD(3);
t65 = t44 * t41 - 0.2e1 * t42 * t62;
t63 = pkin(5) * qJDD(2);
t22 = -t50 * pkin(2) + qJDD(2) * qJ(3) + t66;
t11 = (t44 * t67 - t22 - t63) * t42 + t65;
t60 = t42 * t41 + (t22 + 0.2e1 * t62) * t44;
t12 = -t40 * t67 + t44 * t63 + t60;
t46 = sin(qJ(4));
t48 = cos(qJ(4));
t54 = -t42 * t46 + t44 * t48;
t25 = t54 * qJD(2);
t55 = t42 * t48 + t44 * t46;
t26 = t55 * qJD(2);
t17 = -t25 * mrSges(5,1) + t26 * mrSges(5,2);
t20 = -t26 * qJD(4) + t54 * qJDD(2);
t24 = qJD(4) * mrSges(5,1) - t26 * mrSges(5,3);
t10 = m(5) * (t46 * t11 + t48 * t12) + t20 * mrSges(5,3) - qJDD(4) * mrSges(5,2) + t25 * t17 - qJD(4) * t24;
t56 = -t44 * mrSges(4,1) + t42 * mrSges(4,2);
t53 = qJDD(2) * mrSges(4,3) + t50 * t56;
t21 = t25 * qJD(4) + t55 * qJDD(2);
t23 = -qJD(4) * mrSges(5,2) + t25 * mrSges(5,3);
t9 = m(5) * (t48 * t11 - t46 * t12) - t21 * mrSges(5,3) + qJDD(4) * mrSges(5,1) - t26 * t17 + qJD(4) * t23;
t6 = m(4) * t65 + t46 * t10 + t48 * t9 + (-m(4) * t22 - t53) * t42;
t7 = m(4) * t60 + t48 * t10 + t53 * t44 - t46 * t9;
t61 = m(3) * t41 + t42 * t7 + t44 * t6;
t59 = m(2) * t41 + t61;
t58 = t49 * t31 - t47 * t32;
t57 = qJDD(3) - t58;
t52 = t20 * mrSges(5,1) + t25 * t23 - m(5) * ((-pkin(3) * t44 - pkin(2)) * qJDD(2) + (-t64 * pkin(5) - qJ(3)) * t50 + t57) - t26 * t24 - t21 * mrSges(5,2);
t51 = m(4) * (-qJDD(2) * pkin(2) - t50 * qJ(3) + t57) - t52;
t8 = m(3) * t58 + (-mrSges(3,2) + t68) * t50 + (mrSges(3,1) - t56) * qJDD(2) - t51;
t3 = m(3) * t66 - t50 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t42 * t6 + t44 * t7;
t2 = m(2) * t32 + t49 * t3 - t47 * t8;
t1 = m(2) * t31 + t47 * t3 + t49 * t8;
t4 = [-m(1) * g(1) - t43 * t1 + t45 * t2, t2, t3, t7, t10; -m(1) * g(2) + t45 * t1 + t43 * t2, t1, t8, t6, t9; -m(1) * g(3) + t59, t59, t61, t56 * qJDD(2) - t50 * t68 + t51, -t52;];
f_new = t4;
