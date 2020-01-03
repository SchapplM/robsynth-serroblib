% Calculate vector of cutting forces with Newton-Euler
% S5RPPRP4
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
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
% Datum: 2019-12-31 17:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPPRP4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP4_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP4_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP4_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP4_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP4_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:52:04
% EndTime: 2019-12-31 17:52:05
% DurationCPUTime: 0.39s
% Computational Cost: add. (2372->110), mult. (4168->132), div. (0->0), fcn. (1508->6), ass. (0->53)
t48 = sin(qJ(4));
t50 = cos(qJ(4));
t67 = qJD(1) * qJD(4);
t63 = t50 * t67;
t28 = -t48 * qJDD(1) - t63;
t29 = -t50 * qJDD(1) + t48 * t67;
t68 = qJD(1) * t50;
t35 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t68;
t52 = qJD(1) ^ 2;
t49 = sin(qJ(1));
t51 = cos(qJ(1));
t59 = -t51 * g(1) - t49 * g(2);
t57 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t59;
t78 = -pkin(1) - pkin(2);
t18 = t78 * t52 + t57;
t64 = t49 * g(1) - t51 * g(2);
t54 = -t52 * qJ(2) + qJDD(2) - t64;
t20 = t78 * qJDD(1) + t54;
t46 = sin(pkin(7));
t47 = cos(pkin(7));
t62 = -t46 * t18 + t47 * t20;
t60 = qJDD(1) * pkin(3) - t62;
t69 = qJD(1) * t48;
t31 = qJD(4) * pkin(4) + qJ(5) * t69;
t34 = -qJD(4) * mrSges(6,2) - mrSges(6,3) * t68;
t43 = t50 ^ 2;
t65 = m(6) * (-t31 * t69 - t29 * pkin(4) + qJDD(5) + (-qJ(5) * t43 - pkin(6)) * t52 + t60) + t34 * t68 + t28 * mrSges(6,2);
t32 = qJD(4) * mrSges(6,1) + mrSges(6,3) * t69;
t70 = -qJD(4) * mrSges(5,1) - mrSges(5,3) * t69 - t32;
t82 = (t50 * t35 + t70 * t48) * qJD(1) - (mrSges(5,1) + mrSges(6,1)) * t29 + m(5) * (-t52 * pkin(6) + t60) + t28 * mrSges(5,2) + t65;
t79 = -m(2) - m(3);
t77 = pkin(4) * t52;
t66 = qJD(1) * qJD(5);
t71 = t47 * t18 + t46 * t20;
t14 = -t52 * pkin(3) - qJDD(1) * pkin(6) + t71;
t44 = g(3) + qJDD(3);
t72 = t50 * t14 + t48 * t44;
t76 = t29 * mrSges(6,3) + m(6) * (t29 * qJ(5) - qJD(4) * t31 - t43 * t77 - 0.2e1 * t50 * t66 + t72);
t74 = mrSges(2,1) + mrSges(3,1);
t26 = (mrSges(6,1) * t50 - mrSges(6,2) * t48) * qJD(1);
t37 = t50 * t44;
t61 = t26 * t69 + qJD(4) * t34 + qJDD(4) * mrSges(6,1) + m(6) * (qJDD(4) * pkin(4) + t37 + (-t28 - t63) * qJ(5) + (t50 * t77 - t14 + 0.2e1 * t66) * t48);
t27 = (mrSges(5,1) * t50 - mrSges(5,2) * t48) * qJD(1);
t6 = m(5) * (-t48 * t14 + t37) + qJDD(4) * mrSges(5,1) + t27 * t69 + qJD(4) * t35 + (-mrSges(5,3) - mrSges(6,3)) * t28 + t61;
t7 = m(5) * t72 + t29 * mrSges(5,3) + (-mrSges(5,2) - mrSges(6,2)) * qJDD(4) + t70 * qJD(4) + (-t26 - t27) * t68 + t76;
t4 = m(4) * t71 - t52 * mrSges(4,1) + qJDD(1) * mrSges(4,2) - t48 * t6 + t50 * t7;
t5 = m(4) * t62 - qJDD(1) * mrSges(4,1) - t52 * mrSges(4,2) - t82;
t58 = t47 * t4 - t46 * t5 + m(3) * (-t52 * pkin(1) + t57) + qJDD(1) * mrSges(3,3);
t56 = -m(3) * (-qJDD(1) * pkin(1) + t54) - t46 * t4 - t47 * t5;
t55 = m(4) * t44 + t48 * t7 + t50 * t6;
t2 = m(2) * t64 + (-mrSges(2,2) + mrSges(3,3)) * t52 + t74 * qJDD(1) + t56;
t1 = m(2) * t59 - qJDD(1) * mrSges(2,2) - t74 * t52 + t58;
t3 = [-m(1) * g(1) + t51 * t1 - t49 * t2, t1, -t52 * mrSges(3,1) + t58, t4, t7, -qJDD(4) * mrSges(6,2) - qJD(4) * t32 - t26 * t68 + t76; -m(1) * g(2) + t49 * t1 + t51 * t2, t2, -m(3) * g(3) - t55, t5, t6, -t28 * mrSges(6,3) + t61; (-m(1) + t79) * g(3) - t55, t79 * g(3) - t55, -qJDD(1) * mrSges(3,1) - t52 * mrSges(3,3) - t56, t55, t82, -t29 * mrSges(6,1) - t32 * t69 + t65;];
f_new = t3;
