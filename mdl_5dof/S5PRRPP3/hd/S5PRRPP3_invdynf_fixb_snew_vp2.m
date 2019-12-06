% Calculate vector of cutting forces with Newton-Euler
% S5PRRPP3
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
% Datum: 2019-12-05 16:14
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRRPP3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP3_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP3_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP3_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP3_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP3_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:12:03
% EndTime: 2019-12-05 16:12:05
% DurationCPUTime: 0.62s
% Computational Cost: add. (4742->129), mult. (9701->164), div. (0->0), fcn. (5690->8), ass. (0->67)
t66 = qJD(2) ^ 2;
t60 = sin(pkin(7));
t83 = cos(pkin(7));
t48 = -g(1) * t83 - t60 * g(2);
t58 = -g(3) + qJDD(1);
t62 = sin(qJ(2));
t64 = cos(qJ(2));
t73 = -t62 * t48 + t64 * t58;
t26 = -qJDD(2) * pkin(2) - t66 * pkin(6) - t73;
t61 = sin(qJ(3));
t63 = cos(qJ(3));
t79 = qJD(2) * qJD(3);
t74 = t63 * t79;
t45 = t61 * qJDD(2) + t74;
t75 = t61 * t79;
t46 = t63 * qJDD(2) - t75;
t81 = qJD(2) * t61;
t49 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t81;
t80 = qJD(2) * t63;
t50 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t80;
t59 = sin(pkin(8));
t82 = cos(pkin(8));
t33 = -qJDD(3) * t82 + t59 * t45;
t37 = -qJD(3) * t82 + t59 * t81;
t18 = (-t45 - t74) * qJ(4) + (-t46 + t75) * pkin(3) + t26;
t43 = (-pkin(3) * t63 - qJ(4) * t61) * qJD(2);
t65 = qJD(3) ^ 2;
t84 = t64 * t48 + t62 * t58;
t27 = -t66 * pkin(2) + qJDD(2) * pkin(6) + t84;
t47 = t60 * g(1) - g(2) * t83;
t86 = t63 * t27 - t61 * t47;
t19 = -t65 * pkin(3) + qJDD(3) * qJ(4) + t43 * t80 + t86;
t95 = -2 * qJD(4);
t76 = t59 * t18 + t82 * t19 + t37 * t95;
t38 = t59 * qJD(3) + t81 * t82;
t32 = mrSges(6,1) * t80 + t38 * mrSges(6,2);
t85 = -mrSges(5,1) * t80 - t38 * mrSges(5,3) - t32;
t22 = t37 * mrSges(6,1) - t38 * mrSges(6,3);
t87 = -t37 * mrSges(5,1) - t38 * mrSges(5,2) - t22;
t21 = t37 * pkin(4) - t38 * qJ(5);
t92 = t63 ^ 2 * t66;
t94 = -2 * qJD(5);
t88 = m(6) * (-pkin(4) * t92 - t46 * qJ(5) - t37 * t21 + t80 * t94 + t76) - t46 * mrSges(6,3);
t89 = -mrSges(5,3) - mrSges(6,2);
t8 = m(5) * t76 + t46 * mrSges(5,2) + t33 * t89 + t87 * t37 + t85 * t80 + t88;
t29 = -t37 * mrSges(6,2) - mrSges(6,3) * t80;
t30 = mrSges(5,2) * t80 - t37 * mrSges(5,3);
t34 = t59 * qJDD(3) + t45 * t82;
t70 = t18 * t82 - t59 * t19;
t93 = m(6) * (-qJ(5) * t92 + t46 * pkin(4) + qJDD(5) + ((2 * qJD(4)) + t21) * t38 - t70);
t9 = m(5) * t70 - t93 + (-mrSges(5,1) - mrSges(6,1)) * t46 + (m(5) * t95 + t87) * t38 + t89 * t34 + (-t29 - t30) * t80;
t97 = m(4) * t26 - t46 * mrSges(4,1) + t45 * mrSges(4,2) + (t49 * t61 - t50 * t63) * qJD(2) + t59 * t8 + t82 * t9;
t24 = t61 * t27;
t69 = -qJDD(3) * pkin(3) - t65 * qJ(4) + t43 * t81 + qJDD(4) + t24;
t77 = m(6) * (t33 * pkin(4) - t34 * qJ(5) + t38 * t94 + (t47 + (-pkin(4) * t38 - qJ(5) * t37) * qJD(2)) * t63 + t69) + t37 * t29 + t33 * mrSges(6,1);
t91 = t63 * t47;
t96 = (mrSges(5,2) - mrSges(6,3)) * t34 + t85 * t38 + m(5) * (t69 + t91) + t33 * mrSges(5,1) + t37 * t30 + t77;
t44 = (-mrSges(4,1) * t63 + mrSges(4,2) * t61) * qJD(2);
t10 = m(4) * (-t24 - t91) - t45 * mrSges(4,3) + qJDD(3) * mrSges(4,1) - t44 * t81 + qJD(3) * t50 - t96;
t7 = m(4) * t86 - qJDD(3) * mrSges(4,2) + t46 * mrSges(4,3) - qJD(3) * t49 + t44 * t80 - t59 * t9 + t8 * t82;
t3 = m(3) * t84 - t66 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t61 * t10 + t63 * t7;
t6 = m(3) * t73 + qJDD(2) * mrSges(3,1) - t66 * mrSges(3,2) - t97;
t78 = m(2) * t58 + t62 * t3 + t64 * t6;
t72 = -t63 * t10 - t61 * t7;
t4 = (m(2) + m(3)) * t47 + t72;
t1 = m(2) * t48 + t64 * t3 - t62 * t6;
t2 = [-m(1) * g(1) + t1 * t83 - t60 * t4, t1, t3, t7, t8, -t33 * mrSges(6,2) - t37 * t22 - t32 * t80 + t88; -m(1) * g(2) + t60 * t1 + t4 * t83, t4, t6, t10, t9, -t34 * mrSges(6,3) - t38 * t32 + t77; -m(1) * g(3) + t78, t78, -m(3) * t47 - t72, t97, t96, t46 * mrSges(6,1) + t34 * mrSges(6,2) + t38 * t22 + t29 * t80 + t93;];
f_new = t2;
