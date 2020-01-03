% Calculate vector of cutting forces with Newton-Euler
% S5RRPRP4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRPRP4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP4_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP4_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP4_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP4_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP4_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:52:38
% EndTime: 2019-12-31 19:52:39
% DurationCPUTime: 0.34s
% Computational Cost: add. (2798->108), mult. (3331->133), div. (0->0), fcn. (1408->6), ass. (0->55)
t42 = qJDD(1) + qJDD(2);
t43 = (qJD(1) + qJD(2));
t47 = sin(qJ(1));
t50 = cos(qJ(1));
t63 = t47 * g(1) - g(2) * t50;
t32 = qJDD(1) * pkin(1) + t63;
t52 = qJD(1) ^ 2;
t59 = -g(1) * t50 - g(2) * t47;
t33 = -pkin(1) * t52 + t59;
t46 = sin(qJ(2));
t49 = cos(qJ(2));
t70 = t46 * t32 + t49 * t33;
t80 = -t42 * qJ(3) - (2 * qJD(3) * t43) - t70;
t79 = -m(3) - m(4);
t78 = -pkin(2) - pkin(7);
t41 = t43 ^ 2;
t60 = t49 * t32 - t46 * t33;
t56 = -t41 * qJ(3) + qJDD(3) - t60;
t15 = t78 * t42 + t56;
t45 = sin(qJ(4));
t48 = cos(qJ(4));
t23 = (pkin(4) * t45 - qJ(5) * t48) * t43;
t51 = qJD(4) ^ 2;
t76 = t45 * g(3);
t77 = m(6) * (-qJDD(4) * pkin(4) - t76 - t51 * qJ(5) + qJDD(5) + (t23 * t43 - t15) * t48);
t75 = t43 * t45;
t74 = t43 * t48;
t73 = (mrSges(3,1) - mrSges(4,2));
t72 = -mrSges(3,2) + mrSges(4,3);
t71 = -mrSges(5,3) - mrSges(6,2);
t68 = qJD(4) * t43;
t67 = -m(2) + t79;
t36 = -(qJD(4) * mrSges(6,1)) + mrSges(6,2) * t74;
t62 = -g(3) * t48 + t45 * t15;
t66 = m(6) * (-pkin(4) * t51 + qJDD(4) * qJ(5) + (2 * qJD(5) * qJD(4)) - t23 * t75 + t62) + qJD(4) * t36 + qJDD(4) * mrSges(6,3);
t26 = t42 * t45 + t48 * t68;
t35 = (qJD(4) * mrSges(5,1)) - mrSges(5,3) * t74;
t24 = (mrSges(6,1) * t45 - mrSges(6,3) * t48) * t43;
t61 = t43 * (-t24 - (mrSges(5,1) * t45 + mrSges(5,2) * t48) * t43);
t6 = m(5) * t62 - qJDD(4) * mrSges(5,2) - qJD(4) * t35 + t71 * t26 + t45 * t61 + t66;
t27 = t42 * t48 - t45 * t68;
t34 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t75;
t37 = -mrSges(6,2) * t75 + qJD(4) * mrSges(6,3);
t7 = m(5) * (t15 * t48 + t76) - t77 + t48 * t61 + t71 * t27 + (mrSges(5,1) + mrSges(6,1)) * qJDD(4) + (t34 + t37) * qJD(4);
t64 = -t45 * t7 + t48 * t6;
t58 = -m(4) * (-pkin(2) * t42 + t56) - t45 * t6 - t48 * t7;
t57 = t78 * t41 - t80;
t55 = -t27 * mrSges(6,3) - t36 * t74 + t37 * t75 + t26 * mrSges(6,1) + m(6) * (t26 * pkin(4) - t27 * qJ(5) + (-0.2e1 * qJD(5) * t48 + (pkin(4) * t48 + qJ(5) * t45) * qJD(4)) * t43 + t57);
t54 = m(5) * t57 + t26 * mrSges(5,1) + t27 * mrSges(5,2) + t34 * t75 + t35 * t74 + t55;
t53 = -m(4) * (t41 * pkin(2) + t80) + t54;
t4 = m(3) * t70 - (t73 * t41) + t72 * t42 + t53;
t3 = m(3) * t60 + t72 * t41 + t73 * t42 + t58;
t2 = m(2) * t59 - t52 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t46 * t3 + t49 * t4;
t1 = m(2) * t63 + qJDD(1) * mrSges(2,1) - t52 * mrSges(2,2) + t49 * t3 + t46 * t4;
t5 = [-m(1) * g(1) - t1 * t47 + t2 * t50, t2, t4, -m(4) * g(3) + t64, t6, -t26 * mrSges(6,2) - t24 * t75 + t66; -m(1) * g(2) + t1 * t50 + t2 * t47, t1, t3, -(t41 * mrSges(4,2)) - t42 * mrSges(4,3) - t53, t7, t55; (-m(1) + t67) * g(3) + t64, t67 * g(3) + t64, t79 * g(3) + t64, t42 * mrSges(4,2) - t41 * mrSges(4,3) - t58, t54, -qJDD(4) * mrSges(6,1) + t27 * mrSges(6,2) - qJD(4) * t37 + t24 * t74 + t77;];
f_new = t5;
