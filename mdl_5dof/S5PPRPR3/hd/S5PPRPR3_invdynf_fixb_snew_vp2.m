% Calculate vector of cutting forces with Newton-Euler
% S5PPRPR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
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
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PPRPR3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR3_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR3_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR3_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR3_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR3_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:04:51
% EndTime: 2019-12-05 15:04:52
% DurationCPUTime: 0.37s
% Computational Cost: add. (3699->70), mult. (5700->100), div. (0->0), fcn. (3816->10), ass. (0->48)
t43 = sin(qJ(5));
t45 = cos(qJ(5));
t56 = qJD(3) * qJD(5);
t29 = t43 * qJDD(3) + t45 * t56;
t30 = t45 * qJDD(3) - t43 * t56;
t58 = qJD(3) * t43;
t33 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t58;
t57 = qJD(3) * t45;
t34 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t57;
t47 = qJD(3) ^ 2;
t39 = sin(pkin(7));
t42 = cos(pkin(7));
t32 = -t42 * g(1) - t39 * g(2);
t36 = -g(3) + qJDD(1);
t38 = sin(pkin(8));
t41 = cos(pkin(8));
t25 = t41 * t32 + t38 * t36;
t52 = t39 * g(1) - t42 * g(2);
t31 = qJDD(2) - t52;
t44 = sin(qJ(3));
t46 = cos(qJ(3));
t53 = -t44 * t25 + t46 * t31;
t19 = qJDD(3) * pkin(3) + t53;
t59 = t46 * t25 + t44 * t31;
t20 = -t47 * pkin(3) + t59;
t37 = sin(pkin(9));
t40 = cos(pkin(9));
t51 = t40 * t19 - t37 * t20;
t61 = (t43 * t33 - t45 * t34) * qJD(3) + m(6) * (-qJDD(3) * pkin(4) - t47 * pkin(6) - t51) - t30 * mrSges(6,1) + t29 * mrSges(6,2);
t60 = t37 * t19 + t40 * t20;
t10 = m(5) * t51 + qJDD(3) * mrSges(5,1) - t47 * mrSges(5,2) - t61;
t16 = -t47 * pkin(4) + qJDD(3) * pkin(6) + t60;
t24 = t38 * t32 - t41 * t36;
t23 = qJDD(4) + t24;
t28 = (-mrSges(6,1) * t45 + mrSges(6,2) * t43) * qJD(3);
t13 = m(6) * (-t43 * t16 + t45 * t23) - t29 * mrSges(6,3) + qJDD(5) * mrSges(6,1) - t28 * t58 + qJD(5) * t34;
t14 = m(6) * (t45 * t16 + t43 * t23) + t30 * mrSges(6,3) - qJDD(5) * mrSges(6,2) + t28 * t57 - qJD(5) * t33;
t7 = m(5) * t60 - t47 * mrSges(5,1) - qJDD(3) * mrSges(5,2) - t43 * t13 + t45 * t14;
t5 = m(4) * t53 + qJDD(3) * mrSges(4,1) - t47 * mrSges(4,2) + t40 * t10 + t37 * t7;
t6 = m(4) * t59 - t47 * mrSges(4,1) - qJDD(3) * mrSges(4,2) - t37 * t10 + t40 * t7;
t4 = m(3) * t25 - t44 * t5 + t46 * t6;
t54 = m(5) * t23 + t45 * t13 + t43 * t14;
t9 = (-m(3) - m(4)) * t24 - t54;
t55 = m(2) * t36 + t38 * t4 + t41 * t9;
t49 = m(3) * t31 + t44 * t6 + t46 * t5;
t3 = m(2) * t52 - t49;
t1 = m(2) * t32 - t38 * t9 + t41 * t4;
t2 = [-m(1) * g(1) + t42 * t1 - t39 * t3, t1, t4, t6, t7, t14; -m(1) * g(2) + t39 * t1 + t42 * t3, t3, t9, t5, t10, t13; -m(1) * g(3) + t55, t55, t49, m(4) * t24 + t54, t54, t61;];
f_new = t2;
