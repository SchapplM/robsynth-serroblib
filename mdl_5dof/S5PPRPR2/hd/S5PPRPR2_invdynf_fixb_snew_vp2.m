% Calculate vector of cutting forces with Newton-Euler
% S5PPRPR2
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PPRPR2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR2_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR2_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR2_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR2_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR2_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:03:02
% EndTime: 2019-12-05 15:03:02
% DurationCPUTime: 0.25s
% Computational Cost: add. (1942->71), mult. (3042->93), div. (0->0), fcn. (1818->8), ass. (0->47)
t62 = -pkin(3) - pkin(6);
t61 = mrSges(4,1) - mrSges(5,2);
t60 = -mrSges(4,2) + mrSges(5,3);
t37 = sin(pkin(7));
t39 = cos(pkin(7));
t31 = -t39 * g(1) - t37 * g(2);
t35 = -g(3) + qJDD(1);
t36 = sin(pkin(8));
t38 = cos(pkin(8));
t20 = -t36 * t31 + t38 * t35;
t21 = t38 * t31 + t36 * t35;
t41 = sin(qJ(3));
t43 = cos(qJ(3));
t59 = t41 * t20 + t43 * t21;
t40 = sin(qJ(5));
t58 = qJD(3) * t40;
t42 = cos(qJ(5));
t57 = qJD(3) * t42;
t56 = qJD(3) * qJD(5);
t44 = qJD(3) ^ 2;
t54 = t43 * t20 - t41 * t21;
t48 = -t44 * qJ(4) + qJDD(4) - t54;
t14 = t62 * qJDD(3) + t48;
t27 = (mrSges(6,1) * t40 + mrSges(6,2) * t42) * qJD(3);
t29 = t42 * qJDD(3) - t40 * t56;
t52 = t37 * g(1) - t39 * g(2);
t30 = qJDD(2) - t52;
t32 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t58;
t10 = m(6) * (t42 * t14 - t40 * t30) - t29 * mrSges(6,3) + qJDD(5) * mrSges(6,1) - t27 * t57 + qJD(5) * t32;
t28 = -t40 * qJDD(3) - t42 * t56;
t33 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t57;
t11 = m(6) * (t40 * t14 + t42 * t30) + t28 * mrSges(6,3) - qJDD(5) * mrSges(6,2) - t27 * t58 - qJD(5) * t33;
t49 = -m(5) * (-qJDD(3) * pkin(3) + t48) - t42 * t10 - t40 * t11;
t6 = m(4) * t54 + t61 * qJDD(3) + t60 * t44 + t49;
t46 = qJDD(3) * qJ(4) + 0.2e1 * qJD(4) * qJD(3) + t59;
t50 = -t28 * mrSges(6,1) + m(6) * (t62 * t44 + t46) + t32 * t58 + t33 * t57 + t29 * mrSges(6,2);
t45 = -m(5) * (t44 * pkin(3) - t46) + t50;
t8 = m(4) * t59 + t60 * qJDD(3) - t61 * t44 + t45;
t4 = m(3) * t20 + t41 * t8 + t43 * t6;
t5 = m(3) * t21 - t41 * t6 + t43 * t8;
t55 = m(2) * t35 + t36 * t5 + t38 * t4;
t53 = -m(5) * t30 + t40 * t10 - t42 * t11;
t51 = -m(4) * t30 + t53;
t47 = m(3) * t30 - t51;
t7 = m(2) * t52 - t47;
t1 = m(2) * t31 - t36 * t4 + t38 * t5;
t2 = [-m(1) * g(1) + t39 * t1 - t37 * t7, t1, t5, t8, -t53, t11; -m(1) * g(2) + t37 * t1 + t39 * t7, t7, t4, t6, -t44 * mrSges(5,2) - qJDD(3) * mrSges(5,3) - t45, t10; -m(1) * g(3) + t55, t55, t47, -t51, qJDD(3) * mrSges(5,2) - t44 * mrSges(5,3) - t49, t50;];
f_new = t2;
