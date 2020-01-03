% Calculate vector of cutting forces with Newton-Euler
% S5RPRRP8
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRRP8_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP8_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP8_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP8_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP8_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP8_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP8_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:09
% EndTime: 2019-12-31 18:47:10
% DurationCPUTime: 0.46s
% Computational Cost: add. (3168->112), mult. (4062->136), div. (0->0), fcn. (1455->6), ass. (0->54)
t40 = -qJD(1) + qJD(3);
t38 = t40 ^ 2;
t39 = -qJDD(1) + qJDD(3);
t51 = qJD(1) ^ 2;
t46 = sin(qJ(1));
t49 = cos(qJ(1));
t58 = -t49 * g(1) - t46 * g(2);
t54 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t58;
t76 = -pkin(1) - pkin(2);
t19 = t76 * t51 + t54;
t61 = t46 * g(1) - t49 * g(2);
t52 = -t51 * qJ(2) + qJDD(2) - t61;
t21 = t76 * qJDD(1) + t52;
t45 = sin(qJ(3));
t48 = cos(qJ(3));
t60 = -t45 * t19 + t48 * t21;
t13 = -t39 * pkin(3) - t38 * pkin(7) - t60;
t44 = sin(qJ(4));
t47 = cos(qJ(4));
t63 = qJD(4) * t40;
t27 = t44 * t39 + t47 * t63;
t28 = t47 * t39 - t44 * t63;
t72 = t40 * t44;
t30 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t72;
t31 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t72;
t71 = t40 * t47;
t33 = mrSges(6,2) * t71 + qJD(4) * mrSges(6,3);
t64 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t71 + t33;
t73 = m(6) * (-t28 * pkin(4) - t27 * qJ(5) + (-0.2e1 * qJD(5) * t44 + (pkin(4) * t44 - qJ(5) * t47) * qJD(4)) * t40 + t13) - t28 * mrSges(6,1);
t80 = (-t64 * t47 + (t30 - t31) * t44) * t40 + m(5) * t13 - t28 * mrSges(5,1) + (mrSges(5,2) - mrSges(6,3)) * t27 + t73;
t77 = -m(3) - m(4);
t66 = t48 * t19 + t45 * t21;
t14 = -t38 * pkin(3) + t39 * pkin(7) + t66;
t24 = (-pkin(4) * t47 - qJ(5) * t44) * t40;
t50 = qJD(4) ^ 2;
t74 = t47 * g(3);
t75 = m(6) * (-qJDD(4) * pkin(4) - t74 - t50 * qJ(5) + qJDD(5) + (t24 * t40 + t14) * t44);
t70 = mrSges(2,1) + mrSges(3,1);
t68 = mrSges(5,3) + mrSges(6,2);
t67 = t44 * g(3) + t47 * t14;
t62 = -m(2) + t77;
t25 = (-mrSges(6,1) * t47 - mrSges(6,3) * t44) * t40;
t59 = m(6) * (-t50 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t24 * t71 + t67) + t25 * t71 + qJD(4) * t31 + qJDD(4) * mrSges(6,3);
t26 = (-mrSges(5,1) * t47 + mrSges(5,2) * t44) * t40;
t6 = m(5) * t67 - qJDD(4) * mrSges(5,2) - qJD(4) * t30 + t26 * t71 + t68 * t28 + t59;
t7 = m(5) * (-t44 * t14 + t74) - t75 + (-t25 - t26) * t72 - t68 * t27 + (mrSges(5,1) + mrSges(6,1)) * qJDD(4) + t64 * qJD(4);
t57 = -t44 * t6 - t47 * t7;
t4 = m(4) * t66 - t38 * mrSges(4,1) - t39 * mrSges(4,2) - t44 * t7 + t47 * t6;
t5 = m(4) * t60 + t39 * mrSges(4,1) - t38 * mrSges(4,2) - t80;
t56 = t48 * t4 - t45 * t5 + m(3) * (-t51 * pkin(1) + t54) + qJDD(1) * mrSges(3,3);
t53 = -m(3) * (-qJDD(1) * pkin(1) + t52) - t45 * t4 - t48 * t5;
t2 = m(2) * t61 + (-mrSges(2,2) + mrSges(3,3)) * t51 + t70 * qJDD(1) + t53;
t1 = m(2) * t58 - qJDD(1) * mrSges(2,2) - t70 * t51 + t56;
t3 = [-m(1) * g(1) + t49 * t1 - t46 * t2, t1, -t51 * mrSges(3,1) + t56, t4, t6, t28 * mrSges(6,2) + t59; -m(1) * g(2) + t46 * t1 + t49 * t2, t2, t77 * g(3) + t57, t5, t7, -t27 * mrSges(6,3) + (-t44 * t31 - t47 * t33) * t40 + t73; (-m(1) + t62) * g(3) + t57, t62 * g(3) + t57, -qJDD(1) * mrSges(3,1) - t51 * mrSges(3,3) - t53, m(4) * g(3) - t57, t80, -qJDD(4) * mrSges(6,1) + t27 * mrSges(6,2) - qJD(4) * t33 + t25 * t72 + t75;];
f_new = t3;
