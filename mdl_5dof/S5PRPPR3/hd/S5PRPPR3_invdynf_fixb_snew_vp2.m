% Calculate vector of cutting forces with Newton-Euler
% S5PRPPR3
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRPPR3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR3_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR3_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR3_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR3_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR3_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:26:21
% EndTime: 2019-12-05 15:26:22
% DurationCPUTime: 0.28s
% Computational Cost: add. (2184->78), mult. (3411->99), div. (0->0), fcn. (1818->8), ass. (0->48)
t67 = -pkin(3) - pkin(6);
t66 = mrSges(4,1) - mrSges(5,2);
t65 = -mrSges(4,2) + mrSges(5,3);
t42 = sin(pkin(7));
t44 = cos(pkin(7));
t33 = -g(1) * t44 - g(2) * t42;
t40 = -g(3) + qJDD(1);
t46 = sin(qJ(2));
t48 = cos(qJ(2));
t57 = -t33 * t46 + t48 * t40;
t20 = qJDD(2) * pkin(2) + t57;
t49 = qJD(2) ^ 2;
t63 = t48 * t33 + t46 * t40;
t21 = -pkin(2) * t49 + t63;
t41 = sin(pkin(8));
t43 = cos(pkin(8));
t64 = t41 * t20 + t43 * t21;
t45 = sin(qJ(5));
t62 = qJD(2) * t45;
t47 = cos(qJ(5));
t61 = qJD(2) * t47;
t60 = qJD(2) * qJD(5);
t58 = t20 * t43 - t41 * t21;
t52 = -qJ(4) * t49 + qJDD(4) - t58;
t14 = t67 * qJDD(2) + t52;
t28 = (mrSges(6,1) * t45 + mrSges(6,2) * t47) * qJD(2);
t30 = qJDD(2) * t47 - t45 * t60;
t32 = g(1) * t42 - t44 * g(2);
t31 = qJDD(3) - t32;
t34 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t62;
t10 = m(6) * (t14 * t47 - t31 * t45) - t30 * mrSges(6,3) + qJDD(5) * mrSges(6,1) - t28 * t61 + qJD(5) * t34;
t29 = -qJDD(2) * t45 - t47 * t60;
t35 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t61;
t11 = m(6) * (t14 * t45 + t31 * t47) + t29 * mrSges(6,3) - qJDD(5) * mrSges(6,2) - t28 * t62 - qJD(5) * t35;
t53 = -m(5) * (-qJDD(2) * pkin(3) + t52) - t47 * t10 - t45 * t11;
t6 = m(4) * t58 + t66 * qJDD(2) + t65 * t49 + t53;
t51 = qJDD(2) * qJ(4) + 0.2e1 * qJD(4) * qJD(2) + t64;
t54 = -t29 * mrSges(6,1) + m(6) * (t67 * t49 + t51) + t34 * t62 + t35 * t61 + t30 * mrSges(6,2);
t50 = -m(5) * (pkin(3) * t49 - t51) + t54;
t8 = m(4) * t64 + t65 * qJDD(2) - t66 * t49 + t50;
t4 = m(3) * t57 + qJDD(2) * mrSges(3,1) - t49 * mrSges(3,2) + t41 * t8 + t43 * t6;
t5 = m(3) * t63 - t49 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t41 * t6 + t43 * t8;
t59 = m(2) * t40 + t48 * t4 + t46 * t5;
t56 = -m(5) * t31 + t10 * t45 - t47 * t11;
t55 = -m(4) * t31 + t56;
t7 = (m(2) + m(3)) * t32 + t55;
t1 = m(2) * t33 - t4 * t46 + t48 * t5;
t2 = [-m(1) * g(1) + t1 * t44 - t42 * t7, t1, t5, t8, -t56, t11; -m(1) * g(2) + t1 * t42 + t44 * t7, t7, t4, t6, -t49 * mrSges(5,2) - qJDD(2) * mrSges(5,3) - t50, t10; -m(1) * g(3) + t59, t59, -m(3) * t32 - t55, -t55, qJDD(2) * mrSges(5,2) - t49 * mrSges(5,3) - t53, t54;];
f_new = t2;
