% Calculate vector of cutting forces with Newton-Euler
% S5RRRRP2
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
% Datum: 2019-12-05 18:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRRP2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP2_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP2_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP2_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP2_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:47:35
% EndTime: 2019-12-05 18:47:38
% DurationCPUTime: 0.81s
% Computational Cost: add. (10464->137), mult. (13296->174), div. (0->0), fcn. (7843->8), ass. (0->68)
t62 = qJDD(1) + qJDD(2);
t67 = sin(qJ(3));
t71 = cos(qJ(3));
t64 = qJD(1) + qJD(2);
t87 = qJD(3) * t64;
t46 = t67 * t62 + t71 * t87;
t47 = t71 * t62 - t67 * t87;
t95 = t64 * t67;
t53 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t95;
t94 = t64 * t71;
t54 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t94;
t60 = t64 ^ 2;
t66 = sin(qJ(4));
t70 = cos(qJ(4));
t43 = (t66 * t71 + t67 * t70) * t64;
t26 = -t43 * qJD(4) - t66 * t46 + t70 * t47;
t42 = (-t66 * t67 + t70 * t71) * t64;
t27 = t42 * qJD(4) + t70 * t46 + t66 * t47;
t63 = qJD(3) + qJD(4);
t36 = -t63 * mrSges(6,2) + t42 * mrSges(6,3);
t37 = -t63 * mrSges(5,2) + t42 * mrSges(5,3);
t40 = t63 * mrSges(5,1) - t43 * mrSges(5,3);
t55 = qJD(3) * pkin(3) - pkin(8) * t95;
t65 = t71 ^ 2;
t69 = sin(qJ(1));
t73 = cos(qJ(1));
t88 = t73 * g(2) + t69 * g(3);
t51 = qJDD(1) * pkin(1) + t88;
t74 = qJD(1) ^ 2;
t82 = t69 * g(2) - t73 * g(3);
t52 = -t74 * pkin(1) + t82;
t68 = sin(qJ(2));
t72 = cos(qJ(2));
t80 = t72 * t51 - t68 * t52;
t77 = -t62 * pkin(2) - t80;
t75 = -t47 * pkin(3) + t55 * t95 + (-pkin(8) * t65 - pkin(7)) * t60 + t77;
t38 = t63 * pkin(4) - t43 * qJ(5);
t39 = t63 * mrSges(6,1) - t43 * mrSges(6,3);
t41 = t42 ^ 2;
t84 = m(6) * (-t26 * pkin(4) - t41 * qJ(5) + t43 * t38 + qJDD(5) + t75) + t27 * mrSges(6,2) + t43 * t39;
t76 = m(5) * t75 + t27 * mrSges(5,2) - (mrSges(5,1) + mrSges(6,1)) * t26 + t43 * t40 - (t37 + t36) * t42 + t84;
t103 = (t67 * t53 - t71 * t54) * t64 - t47 * mrSges(4,1) + t46 * mrSges(4,2) + m(4) * (-t60 * pkin(7) + t77) + t76;
t101 = -m(2) - m(3);
t32 = -t42 * mrSges(5,1) + t43 * mrSges(5,2);
t61 = qJDD(3) + qJDD(4);
t31 = -t42 * mrSges(6,1) + t43 * mrSges(6,2);
t89 = t68 * t51 + t72 * t52;
t34 = -t60 * pkin(2) + t62 * pkin(7) + t89;
t93 = t67 * t34;
t98 = pkin(3) * t60;
t18 = qJDD(3) * pkin(3) - t46 * pkin(8) - t93 + (pkin(8) * t87 + t67 * t98 - g(1)) * t71;
t83 = -t67 * g(1) + t71 * t34;
t19 = t47 * pkin(8) - qJD(3) * t55 - t65 * t98 + t83;
t91 = t66 * t18 + t70 * t19;
t85 = m(6) * (-t41 * pkin(4) + t26 * qJ(5) + 0.2e1 * qJD(5) * t42 - t63 * t38 + t91) + t42 * t31 + t26 * mrSges(6,3);
t10 = m(5) * t91 + t26 * mrSges(5,3) + t42 * t32 + (-t40 - t39) * t63 + (-mrSges(5,2) - mrSges(6,2)) * t61 + t85;
t45 = (-mrSges(4,1) * t71 + mrSges(4,2) * t67) * t64;
t81 = t70 * t18 - t66 * t19;
t86 = m(6) * (-0.2e1 * qJD(5) * t43 + (t42 * t63 - t27) * qJ(5) + (t42 * t43 + t61) * pkin(4) + t81) + t63 * t36 + t61 * mrSges(6,1);
t9 = m(5) * t81 + t61 * mrSges(5,1) + t63 * t37 + (-t32 - t31) * t43 + (-mrSges(5,3) - mrSges(6,3)) * t27 + t86;
t6 = m(4) * (-t71 * g(1) - t93) - t46 * mrSges(4,3) + qJDD(3) * mrSges(4,1) - t45 * t95 + qJD(3) * t54 + t66 * t10 + t70 * t9;
t7 = m(4) * t83 - qJDD(3) * mrSges(4,2) + t47 * mrSges(4,3) - qJD(3) * t53 + t70 * t10 + t45 * t94 - t66 * t9;
t100 = t71 * t6 + t67 * t7;
t8 = m(3) * t80 + t62 * mrSges(3,1) - t60 * mrSges(3,2) - t103;
t3 = m(3) * t89 - t60 * mrSges(3,1) - t62 * mrSges(3,2) - t67 * t6 + t71 * t7;
t2 = m(2) * t88 + qJDD(1) * mrSges(2,1) - t74 * mrSges(2,2) + t68 * t3 + t72 * t8;
t1 = m(2) * t82 - t74 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t72 * t3 - t68 * t8;
t4 = [(-m(1) + t101) * g(1) + t100, t1, t3, t7, t10, -t61 * mrSges(6,2) - t63 * t39 + t85; -m(1) * g(2) - t69 * t1 - t73 * t2, t2, t8, t6, t9, -t27 * mrSges(6,3) - t43 * t31 + t86; -m(1) * g(3) + t73 * t1 - t69 * t2, t101 * g(1) + t100, -m(3) * g(1) + t100, t103, t76, -t26 * mrSges(6,1) - t42 * t36 + t84;];
f_new = t4;
