% Calculate vector of cutting forces with Newton-Euler
% S5PPRRP1
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
% Datum: 2019-12-05 15:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PPRRP1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP1_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP1_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP1_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP1_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP1_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:06:41
% EndTime: 2019-12-05 15:06:42
% DurationCPUTime: 0.39s
% Computational Cost: add. (2748->95), mult. (4845->124), div. (0->0), fcn. (2789->8), ass. (0->53)
t52 = sin(qJ(4));
t54 = cos(qJ(4));
t69 = qJD(3) * qJD(4);
t62 = t54 * t69;
t34 = t52 * qJDD(3) + t62;
t35 = t54 * qJDD(3) - t52 * t69;
t71 = qJD(3) * t52;
t41 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t71;
t70 = qJD(3) * t54;
t42 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t70;
t43 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t70;
t56 = qJD(3) ^ 2;
t49 = sin(pkin(7));
t51 = cos(pkin(7));
t38 = -t51 * g(1) - t49 * g(2);
t47 = -g(3) + qJDD(1);
t48 = sin(pkin(8));
t50 = cos(pkin(8));
t23 = -t48 * t38 + t50 * t47;
t24 = t50 * t38 + t48 * t47;
t53 = sin(qJ(3));
t55 = cos(qJ(3));
t61 = t55 * t23 - t53 * t24;
t58 = -qJDD(3) * pkin(3) - t61;
t39 = qJD(4) * pkin(4) - qJ(5) * t71;
t40 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t71;
t46 = t54 ^ 2;
t63 = m(6) * (t39 * t71 - t35 * pkin(4) + qJDD(5) + (-qJ(5) * t46 - pkin(6)) * t56 + t58) + t40 * t71 + t34 * mrSges(6,2);
t80 = -(-t52 * t41 + (t42 + t43) * t54) * qJD(3) - (mrSges(5,1) + mrSges(6,1)) * t35 + m(5) * (-t56 * pkin(6) + t58) + t34 * mrSges(5,2) + t63;
t77 = pkin(4) * t56;
t73 = t53 * t23 + t55 * t24;
t19 = -t56 * pkin(3) + qJDD(3) * pkin(6) + t73;
t60 = t49 * g(1) - t51 * g(2);
t37 = qJDD(2) - t60;
t74 = t54 * t19 + t52 * t37;
t68 = qJD(3) * qJD(5);
t27 = t54 * t37;
t32 = (-mrSges(6,1) * t54 + mrSges(6,2) * t52) * qJD(3);
t33 = (-mrSges(5,1) * t54 + mrSges(5,2) * t52) * qJD(3);
t65 = m(6) * (qJDD(4) * pkin(4) + t27 + (-t34 + t62) * qJ(5) + (t54 * t77 - t19 - 0.2e1 * t68) * t52) + qJD(4) * t42 + qJDD(4) * mrSges(6,1);
t11 = m(5) * (-t52 * t19 + t27) + qJDD(4) * mrSges(5,1) + qJD(4) * t43 + (-mrSges(5,3) - mrSges(6,3)) * t34 + (-t32 - t33) * t71 + t65;
t64 = m(6) * (t35 * qJ(5) - qJD(4) * t39 - t46 * t77 + 0.2e1 * t54 * t68 + t74) + t32 * t70 + t35 * mrSges(6,3);
t12 = m(5) * t74 + t35 * mrSges(5,3) + t33 * t70 + (-mrSges(5,2) - mrSges(6,2)) * qJDD(4) + (-t41 - t40) * qJD(4) + t64;
t6 = m(4) * t73 - t56 * mrSges(4,1) - qJDD(3) * mrSges(4,2) - t52 * t11 + t54 * t12;
t8 = m(4) * t61 + qJDD(3) * mrSges(4,1) - t56 * mrSges(4,2) - t80;
t4 = m(3) * t23 + t53 * t6 + t55 * t8;
t5 = m(3) * t24 - t53 * t8 + t55 * t6;
t67 = m(2) * t47 + t50 * t4 + t48 * t5;
t66 = m(4) * t37 + t54 * t11 + t52 * t12;
t59 = m(3) * t37 + t66;
t7 = m(2) * t60 - t59;
t1 = m(2) * t38 - t48 * t4 + t50 * t5;
t2 = [-m(1) * g(1) + t51 * t1 - t49 * t7, t1, t5, t6, t12, -qJDD(4) * mrSges(6,2) - qJD(4) * t40 + t64; -m(1) * g(2) + t49 * t1 + t51 * t7, t7, t4, t8, t11, -t34 * mrSges(6,3) - t32 * t71 + t65; -m(1) * g(3) + t67, t67, t59, t66, t80, -t35 * mrSges(6,1) - t42 * t70 + t63;];
f_new = t2;
