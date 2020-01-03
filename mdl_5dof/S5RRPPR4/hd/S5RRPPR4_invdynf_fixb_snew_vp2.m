% Calculate vector of cutting forces with Newton-Euler
% S5RRPPR4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRPPR4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR4_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR4_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR4_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR4_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR4_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:27:42
% EndTime: 2019-12-31 19:27:43
% DurationCPUTime: 0.38s
% Computational Cost: add. (4432->87), mult. (5135->109), div. (0->0), fcn. (2168->8), ass. (0->51)
t34 = qJDD(1) + qJDD(2);
t40 = sin(qJ(5));
t43 = cos(qJ(5));
t35 = qJD(1) + qJD(2);
t59 = qJD(5) * t35;
t19 = -t40 * t34 - t43 * t59;
t20 = -t43 * t34 + t40 * t59;
t64 = t35 * t40;
t26 = qJD(5) * mrSges(6,1) + mrSges(6,3) * t64;
t63 = t35 * t43;
t27 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t63;
t42 = sin(qJ(1));
t45 = cos(qJ(1));
t57 = t42 * g(1) - t45 * g(2);
t24 = qJDD(1) * pkin(1) + t57;
t46 = qJD(1) ^ 2;
t54 = -t45 * g(1) - t42 * g(2);
t25 = -t46 * pkin(1) + t54;
t41 = sin(qJ(2));
t44 = cos(qJ(2));
t60 = t41 * t24 + t44 * t25;
t55 = t34 * qJ(3) + 0.2e1 * qJD(3) * t35 + t60;
t65 = -pkin(2) - pkin(3);
t67 = t35 ^ 2;
t13 = t65 * t67 + t55;
t56 = t44 * t24 - t41 * t25;
t47 = -qJ(3) * t67 + qJDD(3) - t56;
t16 = t65 * t34 + t47;
t38 = sin(pkin(8));
t39 = cos(pkin(8));
t53 = -t38 * t13 + t39 * t16;
t68 = -(t40 * t26 - t43 * t27) * t35 + m(6) * (t34 * pkin(4) - pkin(7) * t67 - t53) - t20 * mrSges(6,1) + t19 * mrSges(6,2);
t66 = -m(3) - m(4);
t62 = mrSges(3,1) + mrSges(4,1);
t61 = t39 * t13 + t38 * t16;
t58 = -m(2) + t66;
t11 = -pkin(4) * t67 - t34 * pkin(7) + t61;
t18 = (mrSges(6,1) * t43 - mrSges(6,2) * t40) * t35;
t37 = g(3) + qJDD(4);
t8 = m(6) * (-t40 * t11 + t43 * t37) - t19 * mrSges(6,3) + qJDD(5) * mrSges(6,1) + t18 * t64 + qJD(5) * t27;
t9 = m(6) * (t43 * t11 + t40 * t37) + t20 * mrSges(6,3) - qJDD(5) * mrSges(6,2) - t18 * t63 - qJD(5) * t26;
t6 = m(5) * t61 - mrSges(5,1) * t67 + t34 * mrSges(5,2) - t40 * t8 + t43 * t9;
t7 = m(5) * t53 - t34 * mrSges(5,1) - mrSges(5,2) * t67 - t68;
t51 = -t38 * t7 + t39 * t6 + m(4) * (-pkin(2) * t67 + t55) + t34 * mrSges(4,3);
t50 = -m(4) * (-t34 * pkin(2) + t47) - t38 * t6 - t39 * t7;
t49 = m(5) * t37 + t40 * t9 + t43 * t8;
t4 = m(3) * t56 + t62 * t34 + (-mrSges(3,2) + mrSges(4,3)) * t67 + t50;
t3 = m(3) * t60 - t34 * mrSges(3,2) - t62 * t67 + t51;
t2 = m(2) * t54 - t46 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t44 * t3 - t41 * t4;
t1 = m(2) * t57 + qJDD(1) * mrSges(2,1) - t46 * mrSges(2,2) + t41 * t3 + t44 * t4;
t5 = [-m(1) * g(1) - t42 * t1 + t45 * t2, t2, t3, -mrSges(4,1) * t67 + t51, t6, t9; -m(1) * g(2) + t45 * t1 + t42 * t2, t1, t4, -m(4) * g(3) - t49, t7, t8; (-m(1) + t58) * g(3) - t49, t58 * g(3) - t49, t66 * g(3) - t49, -t34 * mrSges(4,1) - mrSges(4,3) * t67 - t50, t49, t68;];
f_new = t5;
