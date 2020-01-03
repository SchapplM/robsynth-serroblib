% Calculate vector of cutting forces with Newton-Euler
% S4PRPR6
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
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4PRPR6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR6_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR6_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR6_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR6_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR6_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:24:21
% EndTime: 2019-12-31 16:24:22
% DurationCPUTime: 0.32s
% Computational Cost: add. (2243->79), mult. (4831->109), div. (0->0), fcn. (3054->8), ass. (0->50)
t49 = qJD(2) ^ 2;
t44 = cos(pkin(7));
t40 = t44 ^ 2;
t42 = sin(pkin(7));
t64 = t42 ^ 2 + t40;
t68 = t64 * mrSges(4,3);
t67 = pkin(3) * t49;
t43 = sin(pkin(6));
t63 = cos(pkin(6));
t32 = t43 * g(1) - t63 * g(2);
t61 = qJD(2) * qJD(3);
t66 = -t44 * t32 - 0.2e1 * t42 * t61;
t33 = -t63 * g(1) - t43 * g(2);
t41 = -g(3) + qJDD(1);
t46 = sin(qJ(2));
t48 = cos(qJ(2));
t65 = t48 * t33 + t46 * t41;
t62 = pkin(5) * qJDD(2);
t22 = -t49 * pkin(2) + qJDD(2) * qJ(3) + t65;
t11 = (t44 * t67 - t22 - t62) * t42 + t66;
t59 = -t42 * t32 + (t22 + 0.2e1 * t61) * t44;
t12 = -t40 * t67 + t44 * t62 + t59;
t45 = sin(qJ(4));
t47 = cos(qJ(4));
t53 = -t42 * t45 + t44 * t47;
t25 = t53 * qJD(2);
t54 = t42 * t47 + t44 * t45;
t26 = t54 * qJD(2);
t16 = -t25 * mrSges(5,1) + t26 * mrSges(5,2);
t18 = -t26 * qJD(4) + t53 * qJDD(2);
t24 = qJD(4) * mrSges(5,1) - t26 * mrSges(5,3);
t10 = m(5) * (t45 * t11 + t47 * t12) + t18 * mrSges(5,3) - qJDD(4) * mrSges(5,2) + t25 * t16 - qJD(4) * t24;
t55 = -t44 * mrSges(4,1) + t42 * mrSges(4,2);
t52 = qJDD(2) * mrSges(4,3) + t49 * t55;
t19 = t25 * qJD(4) + t54 * qJDD(2);
t23 = -qJD(4) * mrSges(5,2) + t25 * mrSges(5,3);
t9 = m(5) * (t47 * t11 - t45 * t12) - t19 * mrSges(5,3) + qJDD(4) * mrSges(5,1) - t26 * t16 + qJD(4) * t23;
t5 = m(4) * t66 + t45 * t10 + t47 * t9 + (-m(4) * t22 - t52) * t42;
t6 = m(4) * t59 + t47 * t10 + t52 * t44 - t45 * t9;
t3 = m(3) * t65 - t49 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t42 * t5 + t44 * t6;
t58 = -t46 * t33 + t48 * t41;
t56 = qJDD(3) - t58;
t51 = t18 * mrSges(5,1) + t25 * t23 - m(5) * ((-pkin(3) * t44 - pkin(2)) * qJDD(2) + (-t64 * pkin(5) - qJ(3)) * t49 + t56) - t26 * t24 - t19 * mrSges(5,2);
t50 = m(4) * (-qJDD(2) * pkin(2) - t49 * qJ(3) + t56) - t51;
t8 = m(3) * t58 + (-mrSges(3,2) + t68) * t49 + (mrSges(3,1) - t55) * qJDD(2) - t50;
t60 = m(2) * t41 + t46 * t3 + t48 * t8;
t57 = -t42 * t6 - t44 * t5;
t4 = (m(2) + m(3)) * t32 + t57;
t1 = m(2) * t33 + t48 * t3 - t46 * t8;
t2 = [-m(1) * g(1) + t63 * t1 - t43 * t4, t1, t3, t6, t10; -m(1) * g(2) + t43 * t1 + t63 * t4, t4, t8, t5, t9; -m(1) * g(3) + t60, t60, -m(3) * t32 - t57, t55 * qJDD(2) - t49 * t68 + t50, -t51;];
f_new = t2;
