% Calculate vector of cutting forces with Newton-Euler
% S5RPRPP1
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
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
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
% Datum: 2019-12-31 18:09
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRPP1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP1_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP1_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP1_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP1_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP1_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:08:41
% EndTime: 2019-12-31 18:08:42
% DurationCPUTime: 0.68s
% Computational Cost: add. (5955->133), mult. (12482->173), div. (0->0), fcn. (7135->8), ass. (0->65)
t65 = sin(qJ(3));
t67 = cos(qJ(3));
t86 = qJD(1) * qJD(3);
t81 = t67 * t86;
t49 = t65 * qJDD(1) + t81;
t50 = t67 * qJDD(1) - t65 * t86;
t88 = qJD(1) * t65;
t52 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t88;
t87 = qJD(1) * t67;
t53 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t87;
t70 = qJD(1) ^ 2;
t62 = sin(pkin(8));
t89 = cos(pkin(8));
t31 = t62 * t49 - t89 * t50;
t32 = t89 * t49 + t62 * t50;
t39 = t62 * t88 - t89 * t87;
t34 = -qJD(3) * mrSges(5,2) - t39 * mrSges(5,3);
t40 = (t62 * t67 + t89 * t65) * qJD(1);
t35 = qJD(3) * mrSges(5,1) - t40 * mrSges(5,3);
t51 = qJD(3) * pkin(3) - qJ(4) * t88;
t60 = t67 ^ 2;
t66 = sin(qJ(1));
t68 = cos(qJ(1));
t82 = t66 * g(1) - t68 * g(2);
t46 = qJDD(1) * pkin(1) + t82;
t78 = -t68 * g(1) - t66 * g(2);
t48 = -t70 * pkin(1) + t78;
t63 = sin(pkin(7));
t64 = cos(pkin(7));
t79 = t64 * t46 - t63 * t48;
t76 = -qJDD(1) * pkin(2) - t79;
t72 = -t50 * pkin(3) + qJDD(4) + t51 * t88 + (-qJ(4) * t60 - pkin(6)) * t70 + t76;
t36 = -qJD(3) * mrSges(6,1) + t40 * mrSges(6,2);
t37 = -t39 * mrSges(6,2) + qJD(3) * mrSges(6,3);
t74 = t32 * mrSges(6,3) + t40 * t36 - m(6) * (-0.2e1 * qJD(5) * t40 + (qJD(3) * t39 - t32) * qJ(5) + (qJD(3) * t40 + t31) * pkin(4) + t72) - t39 * t37 - t31 * mrSges(6,1);
t73 = m(5) * t72 + t31 * mrSges(5,1) + t32 * mrSges(5,2) + t39 * t34 + t40 * t35 - t74;
t96 = (t65 * t52 - t67 * t53) * qJD(1) + m(4) * (-t70 * pkin(6) + t76) - t50 * mrSges(4,1) + t49 * mrSges(4,2) + t73;
t95 = -2 * qJD(4);
t26 = t39 * pkin(4) - t40 * qJ(5);
t69 = qJD(3) ^ 2;
t90 = t63 * t46 + t64 * t48;
t23 = -t70 * pkin(2) + qJDD(1) * pkin(6) + t90;
t61 = -g(3) + qJDD(2);
t80 = -t65 * t23 + t67 * t61;
t17 = (-t49 + t81) * qJ(4) + (t65 * t67 * t70 + qJDD(3)) * pkin(3) + t80;
t92 = t67 * t23 + t65 * t61;
t18 = -t60 * t70 * pkin(3) + t50 * qJ(4) - qJD(3) * t51 + t92;
t75 = t89 * t17 - t62 * t18;
t94 = m(6) * (-qJDD(3) * pkin(4) - t69 * qJ(5) + qJDD(5) + ((2 * qJD(4)) + t26) * t40 - t75);
t93 = -mrSges(5,3) - mrSges(6,2);
t27 = t39 * mrSges(6,1) - t40 * mrSges(6,3);
t91 = -t39 * mrSges(5,1) - t40 * mrSges(5,2) - t27;
t10 = m(5) * t75 - t94 + (m(5) * t95 + t91) * t40 + t93 * t32 + (mrSges(5,1) + mrSges(6,1)) * qJDD(3) + (t34 + t37) * qJD(3);
t47 = (-mrSges(4,1) * t67 + mrSges(4,2) * t65) * qJD(1);
t83 = t62 * t17 + t89 * t18 + t39 * t95;
t84 = m(6) * (-t69 * pkin(4) + qJDD(3) * qJ(5) + 0.2e1 * qJD(5) * qJD(3) - t39 * t26 + t83) + qJD(3) * t36 + qJDD(3) * mrSges(6,3);
t9 = m(5) * t83 - qJDD(3) * mrSges(5,2) - qJD(3) * t35 + t93 * t31 + t91 * t39 + t84;
t6 = m(4) * t80 + qJDD(3) * mrSges(4,1) - t49 * mrSges(4,3) + qJD(3) * t53 + t89 * t10 - t47 * t88 + t62 * t9;
t7 = m(4) * t92 - qJDD(3) * mrSges(4,2) + t50 * mrSges(4,3) - qJD(3) * t52 - t62 * t10 + t47 * t87 + t89 * t9;
t85 = m(3) * t61 + t67 * t6 + t65 * t7;
t8 = m(3) * t79 + qJDD(1) * mrSges(3,1) - t70 * mrSges(3,2) - t96;
t3 = m(3) * t90 - t70 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t65 * t6 + t67 * t7;
t2 = m(2) * t78 - t70 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t64 * t3 - t63 * t8;
t1 = m(2) * t82 + qJDD(1) * mrSges(2,1) - t70 * mrSges(2,2) + t63 * t3 + t64 * t8;
t4 = [-m(1) * g(1) - t66 * t1 + t68 * t2, t2, t3, t7, t9, -t31 * mrSges(6,2) - t39 * t27 + t84; -m(1) * g(2) + t68 * t1 + t66 * t2, t1, t8, t6, t10, -t74; (-m(1) - m(2)) * g(3) + t85, -m(2) * g(3) + t85, t85, t96, t73, -qJDD(3) * mrSges(6,1) + t32 * mrSges(6,2) - qJD(3) * t37 + t40 * t27 + t94;];
f_new = t4;
