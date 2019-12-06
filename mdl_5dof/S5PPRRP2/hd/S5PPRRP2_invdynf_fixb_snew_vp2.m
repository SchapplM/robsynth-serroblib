% Calculate vector of cutting forces with Newton-Euler
% S5PPRRP2
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
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-05 15:09
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PPRRP2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP2_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP2_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP2_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP2_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP2_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:08:35
% EndTime: 2019-12-05 15:08:36
% DurationCPUTime: 0.40s
% Computational Cost: add. (2723->94), mult. (4705->127), div. (0->0), fcn. (2719->8), ass. (0->51)
t54 = qJD(3) ^ 2;
t46 = sin(pkin(7));
t48 = cos(pkin(7));
t36 = -t48 * g(1) - t46 * g(2);
t44 = -g(3) + qJDD(1);
t45 = sin(pkin(8));
t47 = cos(pkin(8));
t23 = -t45 * t36 + t47 * t44;
t24 = t47 * t36 + t45 * t44;
t50 = sin(qJ(3));
t52 = cos(qJ(3));
t59 = t52 * t23 - t50 * t24;
t18 = -qJDD(3) * pkin(3) - t54 * pkin(6) - t59;
t49 = sin(qJ(4));
t51 = cos(qJ(4));
t62 = qJD(3) * qJD(4);
t32 = t49 * qJDD(3) + t51 * t62;
t33 = t51 * qJDD(3) - t49 * t62;
t64 = qJD(3) * t49;
t37 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t64;
t38 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t64;
t63 = qJD(3) * t51;
t40 = mrSges(6,2) * t63 + qJD(4) * mrSges(6,3);
t65 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t63 + t40;
t69 = m(6) * (-t33 * pkin(4) - t32 * qJ(5) + (-0.2e1 * qJD(5) * t49 + (pkin(4) * t49 - qJ(5) * t51) * qJD(4)) * qJD(3) + t18) - t33 * mrSges(6,1);
t76 = (-t65 * t51 + (t37 - t38) * t49) * qJD(3) + m(5) * t18 - t33 * mrSges(5,1) + (mrSges(5,2) - mrSges(6,3)) * t32 + t69;
t67 = t50 * t23 + t52 * t24;
t19 = -t54 * pkin(3) + qJDD(3) * pkin(6) + t67;
t29 = (-pkin(4) * t51 - qJ(5) * t49) * qJD(3);
t53 = qJD(4) ^ 2;
t57 = t46 * g(1) - t48 * g(2);
t35 = qJDD(2) - t57;
t72 = t51 * t35;
t73 = m(6) * (-qJDD(4) * pkin(4) - t53 * qJ(5) - t72 + qJDD(5) + (qJD(3) * t29 + t19) * t49);
t70 = mrSges(5,3) + mrSges(6,2);
t68 = t51 * t19 + t49 * t35;
t31 = (-mrSges(5,1) * t51 + mrSges(5,2) * t49) * qJD(3);
t30 = (-mrSges(6,1) * t51 - mrSges(6,3) * t49) * qJD(3);
t58 = m(6) * (-t53 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t29 * t63 + t68) + t30 * t63 + qJD(4) * t38 + qJDD(4) * mrSges(6,3);
t11 = m(5) * t68 - qJDD(4) * mrSges(5,2) - qJD(4) * t37 + t31 * t63 + t70 * t33 + t58;
t12 = m(5) * (-t49 * t19 + t72) - t73 - t70 * t32 + (mrSges(5,1) + mrSges(6,1)) * qJDD(4) + t65 * qJD(4) + (-t30 - t31) * t64;
t6 = m(4) * t67 - t54 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t51 * t11 - t49 * t12;
t8 = m(4) * t59 + qJDD(3) * mrSges(4,1) - t54 * mrSges(4,2) - t76;
t4 = m(3) * t23 + t50 * t6 + t52 * t8;
t5 = m(3) * t24 - t50 * t8 + t52 * t6;
t61 = m(2) * t44 + t47 * t4 + t45 * t5;
t60 = m(4) * t35 + t49 * t11 + t51 * t12;
t56 = m(3) * t35 + t60;
t7 = m(2) * t57 - t56;
t1 = m(2) * t36 - t45 * t4 + t47 * t5;
t2 = [-m(1) * g(1) + t48 * t1 - t46 * t7, t1, t5, t6, t11, t33 * mrSges(6,2) + t58; -m(1) * g(2) + t46 * t1 + t48 * t7, t7, t4, t8, t12, -t32 * mrSges(6,3) + (-t49 * t38 - t51 * t40) * qJD(3) + t69; -m(1) * g(3) + t61, t61, t56, t60, t76, -qJDD(4) * mrSges(6,1) + t32 * mrSges(6,2) - qJD(4) * t40 + t30 * t64 + t73;];
f_new = t2;
