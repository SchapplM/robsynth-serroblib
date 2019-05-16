% Calculate vector of cutting forces with Newton-Euler
% S6PPRPRR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
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
% Datum: 2019-05-04 19:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PPRPRR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRPRR1_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRPRR1_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRPRR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_invdynf_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRPRR1_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRPRR1_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRPRR1_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:52:19
% EndTime: 2019-05-04 19:52:24
% DurationCPUTime: 3.18s
% Computational Cost: add. (49084->119), mult. (84739->168), div. (0->0), fcn. (69643->16), ass. (0->81)
t64 = sin(pkin(11));
t69 = cos(pkin(11));
t56 = -t69 * g(1) - t64 * g(2);
t63 = sin(pkin(12));
t68 = cos(pkin(12));
t55 = t64 * g(1) - t69 * g(2);
t61 = -g(3) + qJDD(1);
t66 = sin(pkin(6));
t71 = cos(pkin(6));
t83 = t55 * t71 + t61 * t66;
t37 = -t63 * t56 + t83 * t68;
t47 = -t66 * t55 + t71 * t61 + qJDD(2);
t65 = sin(pkin(7));
t70 = cos(pkin(7));
t105 = t37 * t70 + t47 * t65;
t86 = -t65 * t37 + t70 * t47;
t73 = sin(qJ(5));
t76 = cos(qJ(5));
t52 = (-pkin(5) * t76 - pkin(10) * t73) * qJD(3);
t78 = qJD(5) ^ 2;
t94 = t76 * qJD(3);
t79 = qJD(3) ^ 2;
t38 = t68 * t56 + t83 * t63;
t74 = sin(qJ(3));
t77 = cos(qJ(3));
t85 = t105 * t77 - t74 * t38;
t29 = qJDD(3) * pkin(3) + t85;
t91 = t105 * t74 + t77 * t38;
t30 = -t79 * pkin(3) + t91;
t62 = sin(pkin(13));
t67 = cos(pkin(13));
t96 = t62 * t29 + t67 * t30;
t25 = -t79 * pkin(4) + qJDD(3) * pkin(9) + t96;
t33 = qJDD(4) + t86;
t97 = t76 * t25 + t73 * t33;
t21 = -t78 * pkin(5) + qJDD(5) * pkin(10) + t52 * t94 + t97;
t87 = t67 * t29 - t62 * t30;
t24 = -qJDD(3) * pkin(4) - t79 * pkin(9) - t87;
t93 = qJD(3) * qJD(5);
t88 = t76 * t93;
t53 = t73 * qJDD(3) + t88;
t89 = t73 * t93;
t54 = t76 * qJDD(3) - t89;
t22 = (-t53 - t88) * pkin(10) + (-t54 + t89) * pkin(5) + t24;
t72 = sin(qJ(6));
t75 = cos(qJ(6));
t95 = qJD(3) * t73;
t49 = t75 * qJD(5) - t72 * t95;
t40 = t49 * qJD(6) + t72 * qJDD(5) + t75 * t53;
t50 = t72 * qJD(5) + t75 * t95;
t41 = -t49 * mrSges(7,1) + t50 * mrSges(7,2);
t59 = qJD(6) - t94;
t45 = -t59 * mrSges(7,2) + t49 * mrSges(7,3);
t48 = qJDD(6) - t54;
t18 = m(7) * (-t72 * t21 + t75 * t22) - t40 * mrSges(7,3) + t48 * mrSges(7,1) - t50 * t41 + t59 * t45;
t39 = -t50 * qJD(6) + t75 * qJDD(5) - t72 * t53;
t46 = t59 * mrSges(7,1) - t50 * mrSges(7,3);
t19 = m(7) * (t75 * t21 + t72 * t22) + t39 * mrSges(7,3) - t48 * mrSges(7,2) + t49 * t41 - t59 * t46;
t51 = (-mrSges(6,1) * t76 + mrSges(6,2) * t73) * qJD(3);
t57 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t95;
t15 = m(6) * t97 - qJDD(5) * mrSges(6,2) + t54 * mrSges(6,3) - qJD(5) * t57 - t72 * t18 + t75 * t19 + t51 * t94;
t58 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t94;
t98 = t76 * t33;
t80 = m(7) * (-qJDD(5) * pkin(5) - t78 * pkin(10) - t98 + (qJD(3) * t52 + t25) * t73) - t39 * mrSges(7,1) + t40 * mrSges(7,2) - t49 * t45 + t50 * t46;
t17 = m(6) * (-t73 * t25 + t98) - t53 * mrSges(6,3) + qJDD(5) * mrSges(6,1) - t51 * t95 + qJD(5) * t58 - t80;
t92 = m(5) * t33 + t73 * t15 + t76 * t17;
t12 = m(4) * t86 + t92;
t11 = m(5) * t96 - t79 * mrSges(5,1) - qJDD(3) * mrSges(5,2) + t76 * t15 - t73 * t17;
t103 = m(6) * t24 - t54 * mrSges(6,1) + t53 * mrSges(6,2) + t75 * t18 + t72 * t19 + (t73 * t57 - t76 * t58) * qJD(3);
t13 = m(5) * t87 + qJDD(3) * mrSges(5,1) - t79 * mrSges(5,2) - t103;
t10 = m(4) * t91 - t79 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t67 * t11 - t62 * t13;
t9 = m(4) * t85 + qJDD(3) * mrSges(4,1) - t79 * mrSges(4,2) + t62 * t11 + t67 * t13;
t84 = t10 * t74 + t77 * t9;
t4 = m(3) * t37 - t65 * t12 + t84 * t70;
t8 = m(3) * t38 + t77 * t10 - t74 * t9;
t104 = t4 * t68 + t63 * t8;
t6 = m(3) * t47 + t70 * t12 + t84 * t65;
t90 = m(2) * t61 + t104 * t66 + t71 * t6;
t2 = m(2) * t56 - t63 * t4 + t68 * t8;
t1 = m(2) * t55 + t104 * t71 - t66 * t6;
t3 = [-m(1) * g(1) - t64 * t1 + t69 * t2, t2, t8, t10, t11, t15, t19; -m(1) * g(2) + t69 * t1 + t64 * t2, t1, t4, t9, t13, t17, t18; -m(1) * g(3) + t90, t90, t6, t12, t92, t103, t80;];
f_new  = t3;
