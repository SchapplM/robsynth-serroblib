% Calculate vector of cutting forces with Newton-Euler
% S5PRRPP1
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
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRRPP1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP1_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP1_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP1_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP1_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP1_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:06:14
% EndTime: 2019-12-05 16:06:17
% DurationCPUTime: 0.69s
% Computational Cost: add. (5443->126), mult. (11711->167), div. (0->0), fcn. (7135->8), ass. (0->64)
t65 = sin(qJ(3));
t67 = cos(qJ(3));
t85 = qJD(2) * qJD(3);
t81 = t67 * t85;
t47 = t65 * qJDD(2) + t81;
t48 = t67 * qJDD(2) - t65 * t85;
t87 = qJD(2) * t65;
t52 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t87;
t86 = qJD(2) * t67;
t53 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t86;
t70 = qJD(2) ^ 2;
t62 = sin(pkin(8));
t88 = cos(pkin(8));
t31 = t62 * t47 - t88 * t48;
t32 = t88 * t47 + t62 * t48;
t39 = t62 * t87 - t88 * t86;
t34 = -qJD(3) * mrSges(5,2) - t39 * mrSges(5,3);
t40 = (t62 * t67 + t88 * t65) * qJD(2);
t35 = qJD(3) * mrSges(5,1) - t40 * mrSges(5,3);
t51 = qJD(3) * pkin(3) - qJ(4) * t87;
t60 = t67 ^ 2;
t63 = sin(pkin(7));
t64 = cos(pkin(7));
t49 = t63 * g(1) - t64 * g(2);
t50 = -t64 * g(1) - t63 * g(2);
t66 = sin(qJ(2));
t68 = cos(qJ(2));
t78 = t68 * t49 - t66 * t50;
t76 = -qJDD(2) * pkin(2) - t78;
t72 = -t48 * pkin(3) + qJDD(4) + t51 * t87 + (-qJ(4) * t60 - pkin(6)) * t70 + t76;
t36 = -qJD(3) * mrSges(6,1) + t40 * mrSges(6,2);
t37 = -t39 * mrSges(6,2) + qJD(3) * mrSges(6,3);
t74 = t32 * mrSges(6,3) + t40 * t36 - m(6) * (-0.2e1 * qJD(5) * t40 + (qJD(3) * t39 - t32) * qJ(5) + (qJD(3) * t40 + t31) * pkin(4) + t72) - t39 * t37 - t31 * mrSges(6,1);
t73 = m(5) * t72 + t31 * mrSges(5,1) + t32 * mrSges(5,2) + t39 * t34 + t40 * t35 - t74;
t95 = (t65 * t52 - t67 * t53) * qJD(2) + m(4) * (-t70 * pkin(6) + t76) - t48 * mrSges(4,1) + t47 * mrSges(4,2) + t73;
t94 = -2 * qJD(4);
t23 = t39 * pkin(4) - t40 * qJ(5);
t69 = qJD(3) ^ 2;
t89 = t66 * t49 + t68 * t50;
t30 = -t70 * pkin(2) + qJDD(2) * pkin(6) + t89;
t61 = -g(3) + qJDD(1);
t79 = -t65 * t30 + t67 * t61;
t17 = (-t47 + t81) * qJ(4) + (t65 * t67 * t70 + qJDD(3)) * pkin(3) + t79;
t90 = t67 * t30 + t65 * t61;
t18 = -t60 * t70 * pkin(3) + t48 * qJ(4) - qJD(3) * t51 + t90;
t75 = t88 * t17 - t62 * t18;
t93 = m(6) * (-qJDD(3) * pkin(4) - t69 * qJ(5) + qJDD(5) + ((2 * qJD(4)) + t23) * t40 - t75);
t92 = -mrSges(5,3) - mrSges(6,2);
t24 = t39 * mrSges(6,1) - t40 * mrSges(6,3);
t91 = -t39 * mrSges(5,1) - t40 * mrSges(5,2) - t24;
t10 = m(5) * t75 - t93 + (m(5) * t94 + t91) * t40 + t92 * t32 + (mrSges(5,1) + mrSges(6,1)) * qJDD(3) + (t34 + t37) * qJD(3);
t46 = (-mrSges(4,1) * t67 + mrSges(4,2) * t65) * qJD(2);
t82 = t62 * t17 + t88 * t18 + t39 * t94;
t83 = m(6) * (-t69 * pkin(4) + qJDD(3) * qJ(5) + 0.2e1 * qJD(5) * qJD(3) - t39 * t23 + t82) + qJD(3) * t36 + qJDD(3) * mrSges(6,3);
t9 = m(5) * t82 - qJDD(3) * mrSges(5,2) - qJD(3) * t35 + t92 * t31 + t91 * t39 + t83;
t6 = m(4) * t79 + qJDD(3) * mrSges(4,1) - t47 * mrSges(4,3) + qJD(3) * t53 + t88 * t10 - t46 * t87 + t62 * t9;
t7 = m(4) * t90 - qJDD(3) * mrSges(4,2) + t48 * mrSges(4,3) - qJD(3) * t52 - t62 * t10 + t46 * t86 + t88 * t9;
t84 = m(3) * t61 + t67 * t6 + t65 * t7;
t80 = m(2) * t61 + t84;
t8 = m(3) * t78 + qJDD(2) * mrSges(3,1) - t70 * mrSges(3,2) - t95;
t3 = m(3) * t89 - t70 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t65 * t6 + t67 * t7;
t2 = m(2) * t50 + t68 * t3 - t66 * t8;
t1 = m(2) * t49 + t66 * t3 + t68 * t8;
t4 = [-m(1) * g(1) - t63 * t1 + t64 * t2, t2, t3, t7, t9, -t31 * mrSges(6,2) - t39 * t24 + t83; -m(1) * g(2) + t64 * t1 + t63 * t2, t1, t8, t6, t10, -t74; -m(1) * g(3) + t80, t80, t84, t95, t73, -qJDD(3) * mrSges(6,1) + t32 * mrSges(6,2) - qJD(3) * t37 + t40 * t24 + t93;];
f_new = t4;
