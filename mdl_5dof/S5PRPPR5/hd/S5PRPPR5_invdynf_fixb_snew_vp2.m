% Calculate vector of cutting forces with Newton-Euler
% S5PRPPR5
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% f_new [3x6]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRPPR5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR5_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR5_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR5_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR5_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR5_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:38:04
% EndTime: 2019-12-31 17:38:05
% DurationCPUTime: 0.29s
% Computational Cost: add. (2423->79), mult. (3790->102), div. (0->0), fcn. (1698->8), ass. (0->47)
t45 = sin(qJ(5));
t47 = cos(qJ(5));
t61 = qJD(2) * qJD(5);
t26 = -t45 * qJDD(2) - t47 * t61;
t27 = -t47 * qJDD(2) + t45 * t61;
t63 = qJD(2) * t45;
t32 = qJD(5) * mrSges(6,1) + mrSges(6,3) * t63;
t62 = qJD(2) * t47;
t33 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t62;
t49 = qJD(2) ^ 2;
t42 = sin(pkin(7));
t44 = cos(pkin(7));
t31 = -t44 * g(1) - t42 * g(2);
t40 = -g(3) + qJDD(1);
t46 = sin(qJ(2));
t48 = cos(qJ(2));
t64 = t48 * t31 + t46 * t40;
t57 = qJDD(2) * qJ(3) + 0.2e1 * qJD(3) * qJD(2) + t64;
t67 = -pkin(2) - pkin(3);
t18 = t67 * t49 + t57;
t58 = -t46 * t31 + t48 * t40;
t50 = -t49 * qJ(3) + qJDD(3) - t58;
t20 = t67 * qJDD(2) + t50;
t41 = sin(pkin(8));
t43 = cos(pkin(8));
t55 = -t41 * t18 + t43 * t20;
t68 = -(t45 * t32 - t47 * t33) * qJD(2) + m(6) * (qJDD(2) * pkin(4) - t49 * pkin(6) - t55) - t27 * mrSges(6,1) + t26 * mrSges(6,2);
t66 = mrSges(3,1) + mrSges(4,1);
t65 = t43 * t18 + t41 * t20;
t15 = -t49 * pkin(4) - qJDD(2) * pkin(6) + t65;
t25 = (mrSges(6,1) * t47 - mrSges(6,2) * t45) * qJD(2);
t30 = t42 * g(1) - t44 * g(2);
t28 = qJDD(4) + t30;
t12 = m(6) * (-t45 * t15 + t47 * t28) - t26 * mrSges(6,3) + qJDD(5) * mrSges(6,1) + t25 * t63 + qJD(5) * t33;
t13 = m(6) * (t47 * t15 + t45 * t28) + t27 * mrSges(6,3) - qJDD(5) * mrSges(6,2) - t25 * t62 - qJD(5) * t32;
t7 = m(5) * t65 - t49 * mrSges(5,1) + qJDD(2) * mrSges(5,2) - t45 * t12 + t47 * t13;
t9 = m(5) * t55 - qJDD(2) * mrSges(5,1) - t49 * mrSges(5,2) - t68;
t53 = -t41 * t9 + t43 * t7 + m(4) * (-t49 * pkin(2) + t57) + qJDD(2) * mrSges(4,3);
t4 = m(3) * t64 - qJDD(2) * mrSges(3,2) - t66 * t49 + t53;
t52 = -m(4) * (-qJDD(2) * pkin(2) + t50) - t41 * t7 - t43 * t9;
t5 = m(3) * t58 + (-mrSges(3,2) + mrSges(4,3)) * t49 + t66 * qJDD(2) + t52;
t60 = m(2) * t40 + t46 * t4 + t48 * t5;
t59 = m(5) * t28 + t47 * t12 + t45 * t13;
t56 = m(4) * t30 + t59;
t8 = (m(2) + m(3)) * t30 + t56;
t1 = m(2) * t31 + t48 * t4 - t46 * t5;
t2 = [-m(1) * g(1) + t44 * t1 - t42 * t8, t1, t4, -t49 * mrSges(4,1) + t53, t7, t13; -m(1) * g(2) + t42 * t1 + t44 * t8, t8, t5, -t56, t9, t12; -m(1) * g(3) + t60, t60, -m(3) * t30 - t56, -qJDD(2) * mrSges(4,1) - t49 * mrSges(4,3) - t52, t59, t68;];
f_new = t2;
