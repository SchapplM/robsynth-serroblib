% Calculate vector of cutting forces with Newton-Euler
% S5PPRRP4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
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
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PPRRP4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP4_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP4_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP4_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP4_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP4_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:34:28
% EndTime: 2019-12-31 17:34:28
% DurationCPUTime: 0.30s
% Computational Cost: add. (1650->91), mult. (3027->116), div. (0->0), fcn. (1508->6), ass. (0->49)
t44 = sin(qJ(4));
t46 = cos(qJ(4));
t63 = qJD(3) * qJD(4);
t58 = t46 * t63;
t24 = t44 * qJDD(3) + t58;
t25 = t46 * qJDD(3) - t44 * t63;
t65 = qJD(3) * t44;
t32 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t65;
t64 = qJD(3) * t46;
t33 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t64;
t34 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t64;
t48 = qJD(3) ^ 2;
t42 = sin(pkin(7));
t43 = cos(pkin(7));
t55 = t42 * g(1) - t43 * g(2);
t27 = qJDD(2) - t55;
t29 = -t43 * g(1) - t42 * g(2);
t45 = sin(qJ(3));
t47 = cos(qJ(3));
t57 = t47 * t27 - t45 * t29;
t54 = -qJDD(3) * pkin(3) - t57;
t30 = qJD(4) * pkin(4) - qJ(5) * t65;
t31 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t65;
t40 = t46 ^ 2;
t59 = m(6) * (t30 * t65 - t25 * pkin(4) + qJDD(5) + (-qJ(5) * t40 - pkin(6)) * t48 + t54) + t31 * t65 + t24 * mrSges(6,2);
t74 = -(-t44 * t32 + (t33 + t34) * t46) * qJD(3) - (mrSges(5,1) + mrSges(6,1)) * t25 + m(5) * (-t48 * pkin(6) + t54) + t24 * mrSges(5,2) + t59;
t71 = pkin(4) * t48;
t67 = t45 * t27 + t47 * t29;
t14 = -t48 * pkin(3) + qJDD(3) * pkin(6) + t67;
t41 = g(3) - qJDD(1);
t68 = t46 * t14 + t44 * t41;
t62 = qJD(3) * qJD(5);
t22 = (-mrSges(6,1) * t46 + mrSges(6,2) * t44) * qJD(3);
t61 = t22 * t64 + t25 * mrSges(6,3) + m(6) * (t25 * qJ(5) - qJD(4) * t30 - t40 * t71 + 0.2e1 * t46 * t62 + t68);
t36 = t46 * t41;
t60 = qJD(4) * t33 + qJDD(4) * mrSges(6,1) + m(6) * (qJDD(4) * pkin(4) + t36 + (-t24 + t58) * qJ(5) + (t46 * t71 - t14 - 0.2e1 * t62) * t44);
t23 = (-mrSges(5,1) * t46 + mrSges(5,2) * t44) * qJD(3);
t6 = m(5) * (-t44 * t14 + t36) + qJDD(4) * mrSges(5,1) + qJD(4) * t34 + (-mrSges(5,3) - mrSges(6,3)) * t24 + (-t22 - t23) * t65 + t60;
t7 = m(5) * t68 + t25 * mrSges(5,3) + t23 * t64 + (-mrSges(5,2) - mrSges(6,2)) * qJDD(4) + (-t32 - t31) * qJD(4) + t61;
t4 = m(4) * t67 - t48 * mrSges(4,1) - qJDD(3) * mrSges(4,2) - t44 * t6 + t46 * t7;
t5 = m(4) * t57 + qJDD(3) * mrSges(4,1) - t48 * mrSges(4,2) - t74;
t56 = m(3) * t29 + t47 * t4 - t45 * t5;
t53 = m(3) * t27 + t45 * t4 + t47 * t5;
t52 = m(4) * t41 + t44 * t7 + t46 * t6;
t50 = -m(3) * t41 - t52;
t49 = -m(2) * t41 + t50;
t2 = m(2) * t29 + t56;
t1 = m(2) * t55 - t53;
t3 = [-m(1) * g(1) - t42 * t1 + t43 * t2, t2, t56, t4, t7, -qJDD(4) * mrSges(6,2) - qJD(4) * t31 + t61; -m(1) * g(2) + t43 * t1 + t42 * t2, t1, t50, t5, t6, -t24 * mrSges(6,3) - t22 * t65 + t60; -m(1) * g(3) + t49, t49, t53, t52, t74, -t25 * mrSges(6,1) - t33 * t64 + t59;];
f_new = t3;
