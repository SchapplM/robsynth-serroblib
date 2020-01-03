% Calculate vector of cutting forces with Newton-Euler
% S5RRRRP3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRRP3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP3_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP3_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP3_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP3_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP3_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:49:10
% EndTime: 2019-12-31 21:49:11
% DurationCPUTime: 0.55s
% Computational Cost: add. (6599->111), mult. (6924->141), div. (0->0), fcn. (3373->8), ass. (0->59)
t47 = qJD(1) + qJD(2);
t43 = qJD(3) + t47;
t41 = t43 ^ 2;
t46 = qJDD(1) + qJDD(2);
t42 = qJDD(3) + t46;
t52 = sin(qJ(1));
t56 = cos(qJ(1));
t65 = t52 * g(1) - t56 * g(2);
t38 = qJDD(1) * pkin(1) + t65;
t58 = qJD(1) ^ 2;
t60 = -t56 * g(1) - t52 * g(2);
t39 = -t58 * pkin(1) + t60;
t51 = sin(qJ(2));
t55 = cos(qJ(2));
t62 = t55 * t38 - t51 * t39;
t21 = t46 * pkin(2) + t62;
t45 = t47 ^ 2;
t70 = t51 * t38 + t55 * t39;
t22 = -t45 * pkin(2) + t70;
t50 = sin(qJ(3));
t54 = cos(qJ(3));
t63 = t54 * t21 - t50 * t22;
t16 = -t42 * pkin(3) - t41 * pkin(8) - t63;
t49 = sin(qJ(4));
t53 = cos(qJ(4));
t67 = qJD(4) * t43;
t28 = t49 * t42 + t53 * t67;
t29 = t53 * t42 - t49 * t67;
t76 = t43 * t49;
t34 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t76;
t35 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t76;
t75 = t43 * t53;
t37 = mrSges(6,2) * t75 + qJD(4) * mrSges(6,3);
t68 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t75 + t37;
t72 = m(6) * (-t29 * pkin(4) - t28 * qJ(5) + (-0.2e1 * qJD(5) * t49 + (pkin(4) * t49 - qJ(5) * t53) * qJD(4)) * t43 + t16) - t29 * mrSges(6,1);
t83 = (-t68 * t53 + (t34 - t35) * t49) * t43 + m(5) * t16 - t29 * mrSges(5,1) + (mrSges(5,2) - mrSges(6,3)) * t28 + t72;
t80 = -m(3) - m(4);
t71 = t50 * t21 + t54 * t22;
t17 = -t41 * pkin(3) + t42 * pkin(8) + t71;
t26 = (-mrSges(6,1) * t53 - mrSges(6,3) * t49) * t43;
t27 = (-mrSges(5,1) * t53 + mrSges(5,2) * t49) * t43;
t73 = mrSges(5,3) + mrSges(6,2);
t77 = t53 * g(3);
t25 = (-pkin(4) * t53 - qJ(5) * t49) * t43;
t57 = qJD(4) ^ 2;
t78 = m(6) * (-qJDD(4) * pkin(4) + t77 - t57 * qJ(5) + qJDD(5) + (t25 * t43 + t17) * t49);
t10 = m(5) * (-t49 * t17 - t77) - t78 + (-t26 - t27) * t76 - t73 * t28 + (mrSges(5,1) + mrSges(6,1)) * qJDD(4) + t68 * qJD(4);
t64 = -t49 * g(3) + t53 * t17;
t61 = m(6) * (-t57 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t25 * t75 + t64) + t26 * t75 + qJD(4) * t35 + qJDD(4) * mrSges(6,3);
t9 = m(5) * t64 - qJDD(4) * mrSges(5,2) - qJD(4) * t34 + t27 * t75 + t73 * t29 + t61;
t79 = t53 * t10 + t49 * t9;
t66 = -m(2) + t80;
t6 = m(4) * t63 + t42 * mrSges(4,1) - t41 * mrSges(4,2) - t83;
t5 = m(4) * t71 - t41 * mrSges(4,1) - t42 * mrSges(4,2) - t49 * t10 + t53 * t9;
t4 = m(3) * t70 - t45 * mrSges(3,1) - t46 * mrSges(3,2) + t54 * t5 - t50 * t6;
t3 = m(3) * t62 + t46 * mrSges(3,1) - t45 * mrSges(3,2) + t50 * t5 + t54 * t6;
t2 = m(2) * t60 - t58 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t51 * t3 + t55 * t4;
t1 = m(2) * t65 + qJDD(1) * mrSges(2,1) - t58 * mrSges(2,2) + t55 * t3 + t51 * t4;
t7 = [-m(1) * g(1) - t52 * t1 + t56 * t2, t2, t4, t5, t9, t29 * mrSges(6,2) + t61; -m(1) * g(2) + t56 * t1 + t52 * t2, t1, t3, t6, t10, -t28 * mrSges(6,3) + (-t49 * t35 - t53 * t37) * t43 + t72; (-m(1) + t66) * g(3) + t79, t66 * g(3) + t79, t80 * g(3) + t79, -m(4) * g(3) + t79, t83, -qJDD(4) * mrSges(6,1) + t28 * mrSges(6,2) - qJD(4) * t37 + t26 * t76 + t78;];
f_new = t7;
