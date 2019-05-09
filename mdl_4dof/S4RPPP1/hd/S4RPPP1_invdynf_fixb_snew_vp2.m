% Calculate vector of cutting forces with Newton-Euler
% S4RPPP1
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% f_new [3x5]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 19:07
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4RPPP1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPP1_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPP1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_invdynf_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPP1_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPP1_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPP1_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:07:14
% EndTime: 2019-05-04 19:07:14
% DurationCPUTime: 0.46s
% Computational Cost: add. (1661->136), mult. (4853->202), div. (0->0), fcn. (2829->6), ass. (0->73)
t101 = -2 * qJD(2);
t100 = -2 * qJD(3);
t54 = cos(pkin(4));
t50 = t54 ^ 2;
t99 = pkin(2) * t50;
t98 = mrSges(5,2) * t54;
t55 = sin(qJ(1));
t56 = cos(qJ(1));
t73 = t55 * g(1) - g(2) * t56;
t52 = sin(pkin(4));
t57 = qJD(1) ^ 2;
t94 = t52 * t57;
t29 = qJDD(1) * pkin(1) + qJ(2) * t94 + t73;
t53 = cos(pkin(6));
t97 = t29 * t53;
t51 = sin(pkin(6));
t96 = t51 * t52;
t95 = t52 * t53;
t93 = t54 * t57;
t92 = mrSges(4,1) + mrSges(3,3);
t91 = mrSges(4,2) - mrSges(5,3);
t90 = -mrSges(5,2) - mrSges(4,3);
t68 = -g(1) * t56 - g(2) * t55;
t82 = qJDD(1) * t52;
t30 = -pkin(1) * t57 + qJ(2) * t82 + t68;
t89 = -g(3) * t95 - t51 * t30;
t83 = qJD(1) * t52;
t26 = (mrSges(4,2) * t53 - mrSges(4,3) * t51) * t83;
t88 = t26 + (-mrSges(3,1) * t53 + mrSges(3,2) * t51) * t83;
t77 = mrSges(5,1) * t96;
t34 = (-mrSges(5,3) * t54 + t77) * qJD(1);
t37 = (mrSges(4,1) * t96 + mrSges(4,2) * t54) * qJD(1);
t87 = -t34 - t37;
t62 = -mrSges(4,1) * t95 - t54 * mrSges(4,3);
t35 = t62 * qJD(1);
t36 = (mrSges(5,1) * t95 + t98) * qJD(1);
t86 = t35 - t36;
t85 = qJ(3) * t57;
t48 = t52 ^ 2;
t84 = qJ(4) * t48;
t81 = qJDD(1) * t53;
t31 = (mrSges(3,1) * t54 - mrSges(3,3) * t96) * qJD(1);
t25 = (-pkin(2) * t53 - qJ(3) * t51) * t83;
t65 = -t54 * t51 * t29 + g(3) * t96 - t53 * t30;
t67 = -mrSges(5,2) * t51 - mrSges(5,3) * t53;
t28 = t67 * t83;
t33 = (pkin(3) * t96 - qJ(4) * t54) * qJD(1);
t49 = t53 ^ 2;
t70 = 0.2e1 * qJD(2) * t83;
t61 = t53 * t70 - t65;
t74 = t52 * t81;
t75 = t53 * t83;
t69 = m(5) * (qJDD(4) + (-t49 * t84 - t99) * t57 + (pkin(3) * t95 + qJ(3) * t54) * qJDD(1) + (t25 * t95 + ((2 * qJD(3)) + t33) * t54) * qJD(1) + t61) + qJDD(1) * t98 + mrSges(5,1) * t74 + t54 * qJD(1) * t34 + t28 * t75;
t79 = qJ(3) * qJDD(1);
t60 = m(4) * (t57 * t99 - t54 * t79 + (t54 * t100 + (t101 - t25) * t95) * qJD(1) + t65) - t69;
t5 = m(3) * t61 + ((-mrSges(3,2) + mrSges(4,3)) * t54 + t92 * t95) * qJDD(1) + ((-t31 + t37) * t54 + t88 * t95) * qJD(1) - t60;
t32 = (-t54 * mrSges(3,2) + mrSges(3,3) * t95) * qJD(1);
t76 = t51 * t83;
t58 = t25 * t76 - t50 * t85 + t51 * t70 + qJDD(3) - t89;
t59 = (-pkin(2) - qJ(4)) * qJDD(1) - 0.2e1 * qJD(1) * qJD(4);
t78 = t28 * t76 + qJDD(1) * t77 + m(5) * ((-t53 * t57 * t84 + pkin(3) * t82) * t51 + ((-pkin(3) * t94 - t29) * t53 + t59) * t54 + t58);
t66 = m(4) * ((-pkin(2) * qJDD(1) - t97) * t54 + t58) + t78;
t6 = m(3) * t89 + (-t92 * qJDD(1) + (m(3) * t101 - t88) * qJD(1)) * t96 + (m(3) * t97 + (mrSges(3,1) - t91) * qJDD(1) + (t32 - t86) * qJD(1)) * t54 - t66;
t71 = -t54 * g(3) + qJDD(2);
t64 = pkin(2) * t93 * t96 + t76 * t100 + t71;
t13 = m(5) * (-pkin(3) * t48 * t49 * t57 + (-t29 + (-qJD(1) * t33 - t79) * t51 + (-t54 * t85 + t59) * t53) * t52 + t64);
t72 = t13 + m(4) * ((-pkin(2) * t81 - t29 + (-qJDD(1) * t51 - t53 * t93) * qJ(3)) * t52 + t64) + t35 * t75 + mrSges(4,2) * t74;
t8 = m(3) * t71 + (-m(3) * t29 + ((-mrSges(3,1) - mrSges(5,3)) * qJDD(1) + (-t32 - t36) * qJD(1)) * t53 + ((mrSges(3,2) + t90) * qJDD(1) + (t31 + t87) * qJD(1)) * t51) * t52 + t72;
t80 = t5 * t96 + t54 * t8 + t6 * t95;
t63 = -mrSges(5,3) * qJDD(1) - qJD(1) * t36;
t2 = m(2) * t68 - t57 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t53 * t5 - t51 * t6;
t1 = m(2) * t73 + qJDD(1) * mrSges(2,1) - t57 * mrSges(2,2) - t52 * t8 + (t5 * t51 + t53 * t6) * t54;
t3 = [-m(1) * g(1) - t1 * t55 + t2 * t56, t2, t5, (t63 * t53 + (t87 * qJD(1) + t90 * qJDD(1)) * t51) * t52 + t72, t13 + (t67 * qJDD(1) + (-t34 * t51 - t36 * t53) * qJD(1)) * t52; -m(1) * g(2) + t1 * t56 + t2 * t55, t1, t6, t62 * qJDD(1) + (-t26 * t95 - t54 * t37) * qJD(1) + t60, t63 * t54 + t78; (-m(1) - m(2)) * g(3) + t80, -m(2) * g(3) + t80, t8, (mrSges(4,1) * qJDD(1) + qJD(1) * t26) * t96 + (t86 * qJD(1) + t91 * qJDD(1)) * t54 + t66, t69;];
f_new  = t3;
