% Calculate vector of cutting forces with Newton-Euler
% S6RPPRRP4
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% f_new [3x7]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 14:56
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPPRRP4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP4_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP4_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP4_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP4_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP4_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:54:00
% EndTime: 2019-05-05 14:54:02
% DurationCPUTime: 0.84s
% Computational Cost: add. (9233->150), mult. (16348->179), div. (0->0), fcn. (7917->8), ass. (0->75)
t100 = sin(qJ(5));
t71 = qJD(1) ^ 2;
t102 = -pkin(1) - pkin(2);
t66 = sin(qJ(1));
t69 = cos(qJ(1));
t81 = -t69 * g(1) - t66 * g(2);
t77 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t81;
t38 = t102 * t71 + t77;
t86 = t66 * g(1) - t69 * g(2);
t74 = -t71 * qJ(2) + qJDD(2) - t86;
t40 = t102 * qJDD(1) + t74;
t63 = sin(pkin(9));
t64 = cos(pkin(9));
t82 = -t63 * t38 + t64 * t40;
t21 = qJDD(1) * pkin(3) - t71 * pkin(7) - t82;
t65 = sin(qJ(4));
t68 = cos(qJ(4));
t89 = qJD(1) * qJD(4);
t84 = t68 * t89;
t48 = -t65 * qJDD(1) - t84;
t85 = t65 * t89;
t49 = -t68 * qJDD(1) + t85;
t91 = qJD(1) * t65;
t50 = (qJD(4) * mrSges(5,1)) + mrSges(5,3) * t91;
t90 = t68 * qJD(1);
t51 = -(qJD(4) * mrSges(5,2)) - mrSges(5,3) * t90;
t67 = cos(qJ(5));
t45 = -t100 * qJD(4) + t67 * t91;
t25 = -t45 * qJD(5) - t67 * qJDD(4) + t100 * t48;
t53 = qJD(5) + t90;
t35 = t53 * mrSges(6,1) + t45 * mrSges(6,3);
t43 = qJDD(5) - t49;
t44 = t67 * qJD(4) + t100 * t91;
t104 = 2 * qJD(6);
t28 = -t44 * pkin(5) + t45 * qJ(6);
t36 = -t53 * mrSges(7,1) - t45 * mrSges(7,2);
t52 = t53 ^ 2;
t15 = (-t48 + t84) * pkin(8) + (-t49 - t85) * pkin(4) + t21;
t47 = (pkin(4) * t68 + pkin(8) * t65) * qJD(1);
t70 = qJD(4) ^ 2;
t93 = t64 * t38 + t63 * t40;
t22 = -t71 * pkin(3) - qJDD(1) * pkin(7) + t93;
t61 = g(3) + qJDD(3);
t95 = t68 * t22 + t65 * t61;
t18 = -t70 * pkin(4) + qJDD(4) * pkin(8) - t47 * t90 + t95;
t96 = t100 * t15 + t67 * t18;
t88 = m(7) * (-t52 * pkin(5) + t43 * qJ(6) + t53 * t104 + t44 * t28 + t96) + t53 * t36 + t43 * mrSges(7,3);
t29 = -t44 * mrSges(7,1) + t45 * mrSges(7,3);
t94 = -t44 * mrSges(6,1) - t45 * mrSges(6,2) + t29;
t97 = -mrSges(6,3) - mrSges(7,2);
t8 = m(6) * t96 - t43 * mrSges(6,2) + t97 * t25 - t53 * t35 + t94 * t44 + t88;
t78 = -t100 * t18 + t67 * t15;
t101 = m(7) * (-t43 * pkin(5) - t52 * qJ(6) - t45 * t28 + qJDD(6) - t78);
t26 = t44 * qJD(5) + t100 * qJDD(4) + t67 * t48;
t34 = -t53 * mrSges(6,2) + t44 * mrSges(6,3);
t37 = t44 * mrSges(7,2) + t53 * mrSges(7,3);
t9 = m(6) * t78 - t101 + (t34 + t37) * t53 + t94 * t45 + (mrSges(6,1) + mrSges(7,1)) * t43 + t97 * t26;
t106 = m(5) * t21 - t49 * mrSges(5,1) + t48 * mrSges(5,2) - (t50 * t65 - t51 * t68) * qJD(1) + t100 * t8 + t67 * t9;
t83 = -t65 * t22 + t68 * t61;
t17 = -qJDD(4) * pkin(4) - t70 * pkin(8) - t47 * t91 - t83;
t87 = t44 * t37 - m(7) * (t45 * t104 + (-t44 * t53 - t26) * qJ(6) + (-t45 * t53 + t25) * pkin(5) + t17) - t25 * mrSges(7,1);
t105 = m(6) * t17 + t25 * mrSges(6,1) + (mrSges(6,2) - mrSges(7,3)) * t26 - t44 * t34 - (t35 - t36) * t45 - t87;
t103 = -m(2) - m(3);
t99 = mrSges(2,1) + mrSges(3,1);
t46 = (mrSges(5,1) * t68 - mrSges(5,2) * t65) * qJD(1);
t6 = m(5) * t95 - qJDD(4) * mrSges(5,2) + t49 * mrSges(5,3) - qJD(4) * t50 - t100 * t9 - t46 * t90 + t67 * t8;
t7 = m(5) * t83 + qJDD(4) * mrSges(5,1) - t48 * mrSges(5,3) + qJD(4) * t51 + t46 * t91 - t105;
t4 = m(4) * t93 - t71 * mrSges(4,1) + qJDD(1) * mrSges(4,2) + t68 * t6 - t65 * t7;
t5 = m(4) * t82 - qJDD(1) * mrSges(4,1) - t71 * mrSges(4,2) - t106;
t79 = t64 * t4 - t63 * t5 + m(3) * (-t71 * pkin(1) + t77) + qJDD(1) * mrSges(3,3);
t76 = -m(3) * (-qJDD(1) * pkin(1) + t74) - t63 * t4 - t64 * t5;
t75 = m(4) * t61 + t65 * t6 + t68 * t7;
t2 = m(2) * t86 + (-mrSges(2,2) + mrSges(3,3)) * t71 + t99 * qJDD(1) + t76;
t1 = m(2) * t81 - qJDD(1) * mrSges(2,2) - t99 * t71 + t79;
t3 = [-m(1) * g(1) + t69 * t1 - t66 * t2, t1, -t71 * mrSges(3,1) + t79, t4, t6, t8, -t25 * mrSges(7,2) + t44 * t29 + t88; -m(1) * g(2) + t66 * t1 + t69 * t2, t2, -m(3) * g(3) - t75, t5, t7, t9, -t26 * mrSges(7,3) + t45 * t36 - t87; (-m(1) + t103) * g(3) - t75, t103 * g(3) - t75, -qJDD(1) * mrSges(3,1) - t71 * mrSges(3,3) - t76, t75, t106, t105, -t43 * mrSges(7,1) + t26 * mrSges(7,2) - t45 * t29 - t53 * t37 + t101;];
f_new  = t3;
