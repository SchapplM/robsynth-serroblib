% Calculate vector of cutting forces with Newton-Euler
% S5RPRRR7
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRRR7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR7_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR7_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR7_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR7_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR7_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:03:17
% EndTime: 2019-12-31 19:03:19
% DurationCPUTime: 1.15s
% Computational Cost: add. (13681->137), mult. (26293->181), div. (0->0), fcn. (16460->10), ass. (0->76)
t77 = qJD(1) ^ 2;
t71 = sin(qJ(1));
t75 = cos(qJ(1));
t87 = t71 * g(1) - t75 * g(2);
t51 = qJDD(1) * pkin(1) + t87;
t82 = -t75 * g(1) - t71 * g(2);
t53 = -t77 * pkin(1) + t82;
t66 = sin(pkin(9));
t67 = cos(pkin(9));
t83 = t67 * t51 - t66 * t53;
t32 = -qJDD(1) * pkin(2) - t77 * pkin(6) - t83;
t70 = sin(qJ(3));
t74 = cos(qJ(3));
t89 = qJD(1) * qJD(3);
t86 = t74 * t89;
t55 = t70 * qJDD(1) + t86;
t62 = t70 * t89;
t56 = t74 * qJDD(1) - t62;
t91 = qJD(1) * t70;
t57 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t91;
t90 = t74 * qJD(1);
t58 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t90;
t69 = sin(qJ(4));
t73 = cos(qJ(4));
t49 = t73 * qJD(3) - t69 * t91;
t35 = t49 * qJD(4) + t69 * qJDD(3) + t73 * t55;
t48 = qJDD(4) - t56;
t50 = t69 * qJD(3) + t73 * t91;
t60 = qJD(4) - t90;
t22 = (-t55 - t86) * pkin(7) + (-t56 + t62) * pkin(3) + t32;
t54 = (-pkin(3) * t74 - pkin(7) * t70) * qJD(1);
t76 = qJD(3) ^ 2;
t92 = t66 * t51 + t67 * t53;
t33 = -t77 * pkin(2) + qJDD(1) * pkin(6) + t92;
t65 = -g(3) + qJDD(2);
t93 = t74 * t33 + t70 * t65;
t26 = -t76 * pkin(3) + qJDD(3) * pkin(7) + t54 * t90 + t93;
t85 = t73 * t22 - t69 * t26;
t13 = (t49 * t60 - t35) * pkin(8) + (t49 * t50 + t48) * pkin(4) + t85;
t34 = -t50 * qJD(4) + t73 * qJDD(3) - t69 * t55;
t41 = t60 * pkin(4) - t50 * pkin(8);
t47 = t49 ^ 2;
t94 = t69 * t22 + t73 * t26;
t14 = -t47 * pkin(4) + t34 * pkin(8) - t60 * t41 + t94;
t68 = sin(qJ(5));
t72 = cos(qJ(5));
t36 = t72 * t49 - t68 * t50;
t19 = t36 * qJD(5) + t68 * t34 + t72 * t35;
t37 = t68 * t49 + t72 * t50;
t27 = -t36 * mrSges(6,1) + t37 * mrSges(6,2);
t59 = qJD(5) + t60;
t28 = -t59 * mrSges(6,2) + t36 * mrSges(6,3);
t46 = qJDD(5) + t48;
t11 = m(6) * (t72 * t13 - t68 * t14) - t19 * mrSges(6,3) + t46 * mrSges(6,1) - t37 * t27 + t59 * t28;
t18 = -t37 * qJD(5) + t72 * t34 - t68 * t35;
t29 = t59 * mrSges(6,1) - t37 * mrSges(6,3);
t12 = m(6) * (t68 * t13 + t72 * t14) + t18 * mrSges(6,3) - t46 * mrSges(6,2) + t36 * t27 - t59 * t29;
t38 = -t49 * mrSges(5,1) + t50 * mrSges(5,2);
t39 = -t60 * mrSges(5,2) + t49 * mrSges(5,3);
t7 = m(5) * t85 + t48 * mrSges(5,1) - t35 * mrSges(5,3) + t72 * t11 + t68 * t12 - t50 * t38 + t60 * t39;
t40 = t60 * mrSges(5,1) - t50 * mrSges(5,3);
t8 = m(5) * t94 - t48 * mrSges(5,2) + t34 * mrSges(5,3) - t68 * t11 + t72 * t12 + t49 * t38 - t60 * t40;
t95 = m(4) * t32 - t56 * mrSges(4,1) + t55 * mrSges(4,2) + t69 * t8 + t73 * t7 + (t70 * t57 - t74 * t58) * qJD(1);
t52 = (-mrSges(4,1) * t74 + mrSges(4,2) * t70) * qJD(1);
t84 = -t70 * t33 + t74 * t65;
t25 = -qJDD(3) * pkin(3) - t76 * pkin(7) + t54 * t91 - t84;
t80 = t18 * mrSges(6,1) + t36 * t28 - m(6) * (-t34 * pkin(4) - t47 * pkin(8) + t50 * t41 + t25) - t19 * mrSges(6,2) - t37 * t29;
t78 = m(5) * t25 - t34 * mrSges(5,1) + t35 * mrSges(5,2) - t49 * t39 + t50 * t40 - t80;
t10 = m(4) * t84 + qJDD(3) * mrSges(4,1) - t55 * mrSges(4,3) + qJD(3) * t58 - t52 * t91 - t78;
t6 = m(4) * t93 - qJDD(3) * mrSges(4,2) + t56 * mrSges(4,3) - qJD(3) * t57 + t52 * t90 - t69 * t7 + t73 * t8;
t88 = m(3) * t65 + t74 * t10 + t70 * t6;
t4 = m(3) * t83 + qJDD(1) * mrSges(3,1) - t77 * mrSges(3,2) - t95;
t3 = m(3) * t92 - t77 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t70 * t10 + t74 * t6;
t2 = m(2) * t82 - t77 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t67 * t3 - t66 * t4;
t1 = m(2) * t87 + qJDD(1) * mrSges(2,1) - t77 * mrSges(2,2) + t66 * t3 + t67 * t4;
t5 = [-m(1) * g(1) - t71 * t1 + t75 * t2, t2, t3, t6, t8, t12; -m(1) * g(2) + t75 * t1 + t71 * t2, t1, t4, t10, t7, t11; (-m(1) - m(2)) * g(3) + t88, -m(2) * g(3) + t88, t88, t95, t78, -t80;];
f_new = t5;
