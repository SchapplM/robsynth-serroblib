% Calculate vector of cutting forces with Newton-Euler
% S5RRPRP9
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 20:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRPRP9_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP9_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP9_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP9_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP9_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP9_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP9_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP9_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP9_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:05:45
% EndTime: 2019-12-31 20:05:47
% DurationCPUTime: 1.16s
% Computational Cost: add. (11788->158), mult. (25326->204), div. (0->0), fcn. (16414->8), ass. (0->77)
t105 = cos(qJ(4));
t78 = sin(qJ(2));
t100 = qJD(1) * t78;
t75 = sin(pkin(8));
t76 = cos(pkin(8));
t61 = t76 * qJD(2) - t75 * t100;
t62 = t75 * qJD(2) + t76 * t100;
t44 = -t61 * mrSges(4,1) + t62 * mrSges(4,2);
t80 = cos(qJ(2));
t99 = t80 * qJD(1);
t46 = mrSges(4,2) * t99 + t61 * mrSges(4,3);
t98 = qJD(1) * qJD(2);
t93 = t80 * t98;
t66 = t78 * qJDD(1) + t93;
t49 = t75 * qJDD(2) + t76 * t66;
t72 = t78 * t98;
t67 = t80 * qJDD(1) - t72;
t77 = sin(qJ(4));
t42 = -t105 * t61 + t77 * t62;
t43 = t105 * t62 + t77 * t61;
t28 = t42 * mrSges(6,1) - t43 * mrSges(6,3);
t102 = -t42 * mrSges(5,1) - t43 * mrSges(5,2) - t28;
t83 = qJD(1) ^ 2;
t79 = sin(qJ(1));
t81 = cos(qJ(1));
t95 = t79 * g(1) - t81 * g(2);
t57 = -qJDD(1) * pkin(1) - t83 * pkin(6) - t95;
t32 = (-t66 - t93) * qJ(3) + (-t67 + t72) * pkin(2) + t57;
t64 = (-pkin(2) * t80 - qJ(3) * t78) * qJD(1);
t82 = qJD(2) ^ 2;
t91 = -t81 * g(1) - t79 * g(2);
t58 = -t83 * pkin(1) + qJDD(1) * pkin(6) + t91;
t94 = -t78 * g(3) + t80 * t58;
t36 = -t82 * pkin(2) + qJDD(2) * qJ(3) + t64 * t99 + t94;
t92 = -0.2e1 * qJD(3) * t62 + t76 * t32 - t75 * t36;
t16 = (-t61 * t99 - t49) * pkin(7) + (t61 * t62 - t67) * pkin(3) + t92;
t48 = t76 * qJDD(2) - t75 * t66;
t50 = -pkin(3) * t99 - t62 * pkin(7);
t60 = t61 ^ 2;
t96 = 0.2e1 * qJD(3) * t61 + t75 * t32 + t76 * t36;
t18 = -t60 * pkin(3) + t48 * pkin(7) + t50 * t99 + t96;
t103 = t105 * t18 + t77 * t16;
t104 = -mrSges(5,3) - mrSges(6,2);
t23 = t43 * qJD(4) - t105 * t48 + t77 * t49;
t71 = qJD(4) - t99;
t38 = t71 * mrSges(5,1) - t43 * mrSges(5,3);
t63 = qJDD(4) - t67;
t27 = t42 * pkin(4) - t43 * qJ(5);
t39 = -t71 * mrSges(6,1) + t43 * mrSges(6,2);
t70 = t71 ^ 2;
t97 = m(6) * (-t70 * pkin(4) + t63 * qJ(5) + 0.2e1 * qJD(5) * t71 - t42 * t27 + t103) + t71 * t39 + t63 * mrSges(6,3);
t7 = m(5) * t103 - t63 * mrSges(5,2) + t102 * t42 + t104 * t23 - t71 * t38 + t97;
t89 = t105 * t16 - t77 * t18;
t106 = m(6) * (-t63 * pkin(4) - t70 * qJ(5) + t43 * t27 + qJDD(5) - t89);
t24 = -t42 * qJD(4) + t105 * t49 + t77 * t48;
t37 = -t71 * mrSges(5,2) - t42 * mrSges(5,3);
t40 = -t42 * mrSges(6,2) + t71 * mrSges(6,3);
t8 = m(5) * t89 - t106 + (t37 + t40) * t71 + (mrSges(5,1) + mrSges(6,1)) * t63 + t102 * t43 + t104 * t24;
t5 = m(4) * t92 - t67 * mrSges(4,1) - t49 * mrSges(4,3) + t105 * t8 - t62 * t44 - t46 * t99 + t77 * t7;
t47 = -mrSges(4,1) * t99 - t62 * mrSges(4,3);
t6 = m(4) * t96 + t67 * mrSges(4,2) + t48 * mrSges(4,3) + t105 * t7 + t61 * t44 + t47 * t99 - t77 * t8;
t68 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t100;
t69 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t99;
t108 = m(3) * t57 - t67 * mrSges(3,1) + t66 * mrSges(3,2) + t76 * t5 + t75 * t6 + (t78 * t68 - t80 * t69) * qJD(1);
t101 = -t80 * g(3) - t78 * t58;
t65 = (-mrSges(3,1) * t80 + mrSges(3,2) * t78) * qJD(1);
t35 = -qJDD(2) * pkin(2) - t82 * qJ(3) + t64 * t100 + qJDD(3) - t101;
t86 = -t48 * pkin(3) - t60 * pkin(7) + t62 * t50 + t35;
t88 = t24 * mrSges(6,3) + t43 * t39 - m(6) * (-0.2e1 * qJD(5) * t43 + (t42 * t71 - t24) * qJ(5) + (t43 * t71 + t23) * pkin(4) + t86) - t23 * mrSges(6,1) - t42 * t40;
t85 = m(5) * t86 + t23 * mrSges(5,1) + t24 * mrSges(5,2) + t42 * t37 + t43 * t38 - t88;
t84 = m(4) * t35 - t48 * mrSges(4,1) + t49 * mrSges(4,2) - t61 * t46 + t62 * t47 + t85;
t10 = m(3) * t101 + qJDD(2) * mrSges(3,1) - t66 * mrSges(3,3) + qJD(2) * t69 - t65 * t100 - t84;
t4 = m(3) * t94 - qJDD(2) * mrSges(3,2) + t67 * mrSges(3,3) - qJD(2) * t68 - t75 * t5 + t76 * t6 + t65 * t99;
t107 = t80 * t10 + t78 * t4;
t2 = m(2) * t95 + qJDD(1) * mrSges(2,1) - t83 * mrSges(2,2) - t108;
t1 = m(2) * t91 - t83 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t78 * t10 + t80 * t4;
t3 = [-m(1) * g(1) + t81 * t1 - t79 * t2, t1, t4, t6, t7, -t23 * mrSges(6,2) - t42 * t28 + t97; -m(1) * g(2) + t79 * t1 + t81 * t2, t2, t10, t5, t8, -t88; (-m(1) - m(2)) * g(3) + t107, -m(2) * g(3) + t107, t108, t84, t85, -t63 * mrSges(6,1) + t24 * mrSges(6,2) + t43 * t28 - t71 * t40 + t106;];
f_new = t3;
