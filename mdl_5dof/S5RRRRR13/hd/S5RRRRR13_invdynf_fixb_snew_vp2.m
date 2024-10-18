% Calculate vector of cutting forces with Newton-Euler
% S5RRRRR13
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
% m [6x1]
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
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRRR13_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR13_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR13_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR13_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR13_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR13_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR13_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 17:30:36
% EndTime: 2024-09-27 17:30:38
% DurationCPUTime: 1.47s
% Computational Cost: add. (35263->123), mult. (41197->176), div. (0->0), fcn. (26143->12), ass. (0->80)
t99 = -m(3) - m(4);
t63 = sin(pkin(5));
t98 = pkin(9) * t63;
t64 = cos(pkin(5));
t97 = t64 * g(3);
t61 = qJD(1) + qJD(2);
t56 = qJD(3) + t61;
t54 = t56 ^ 2;
t59 = qJDD(1) + qJDD(2);
t55 = qJDD(3) + t59;
t69 = sin(qJ(1));
t74 = cos(qJ(1));
t82 = t69 * g(1) - t74 * g(2);
t47 = qJDD(1) * pkin(1) + t82;
t75 = qJD(1) ^ 2;
t78 = -t74 * g(1) - t69 * g(2);
t48 = -t75 * pkin(1) + t78;
t68 = sin(qJ(2));
t73 = cos(qJ(2));
t79 = t73 * t47 - t68 * t48;
t33 = t59 * pkin(2) + t79;
t58 = t61 ^ 2;
t89 = t68 * t47 + t73 * t48;
t34 = -t58 * pkin(2) + t89;
t67 = sin(qJ(3));
t72 = cos(qJ(3));
t80 = t72 * t33 - t67 * t34;
t21 = t55 * pkin(3) + t54 * t98 + t80;
t96 = t21 * t64;
t95 = t54 * t63 ^ 2;
t94 = t56 * t63;
t66 = sin(qJ(4));
t93 = t63 * t66;
t71 = cos(qJ(4));
t92 = t63 * t71;
t88 = qJD(4) * t56;
t41 = (t55 * t66 + t71 * t88) * t63;
t50 = t64 * t55 + qJDD(4);
t51 = t64 * t56 + qJD(4);
t90 = t67 * t33 + t72 * t34;
t22 = -t54 * pkin(3) + t55 * t98 + t90;
t81 = -t66 * t22 + t71 * t96;
t15 = t50 * pkin(4) - t41 * pkin(10) + (pkin(4) * t66 * t95 + (pkin(10) * t51 * t56 - g(3)) * t63) * t71 + t81;
t84 = t56 * t93;
t40 = t51 * pkin(4) - pkin(10) * t84;
t42 = (t55 * t71 - t66 * t88) * t63;
t77 = -g(3) * t93 + t71 * t22 + t66 * t96;
t85 = t71 ^ 2 * t95;
t16 = -pkin(4) * t85 + t42 * pkin(10) - t51 * t40 + t77;
t65 = sin(qJ(5));
t70 = cos(qJ(5));
t35 = (-t65 * t66 + t70 * t71) * t94;
t25 = t35 * qJD(5) + t70 * t41 + t65 * t42;
t36 = (t65 * t71 + t66 * t70) * t94;
t27 = -t35 * mrSges(6,1) + t36 * mrSges(6,2);
t49 = qJD(5) + t51;
t28 = -t49 * mrSges(6,2) + t35 * mrSges(6,3);
t46 = qJDD(5) + t50;
t11 = m(6) * (t70 * t15 - t65 * t16) - t25 * mrSges(6,3) + t46 * mrSges(6,1) - t36 * t27 + t49 * t28;
t24 = -t36 * qJD(5) - t65 * t41 + t70 * t42;
t29 = t49 * mrSges(6,1) - t36 * mrSges(6,3);
t12 = m(6) * (t65 * t15 + t70 * t16) + t24 * mrSges(6,3) - t46 * mrSges(6,2) + t35 * t27 - t49 * t29;
t37 = t51 * mrSges(5,1) - mrSges(5,3) * t84;
t39 = (-mrSges(5,1) * t71 + mrSges(5,2) * t66) * t94;
t83 = t56 * t92;
t10 = m(5) * t77 - t50 * mrSges(5,2) + t42 * mrSges(5,3) - t65 * t11 + t70 * t12 - t51 * t37 + t39 * t83;
t91 = t66 * t10;
t87 = -m(2) + t99;
t38 = -t51 * mrSges(5,2) + mrSges(5,3) * t83;
t76 = -t24 * mrSges(6,1) - t35 * t28 + m(6) * (-pkin(10) * t85 - t42 * pkin(4) - t97 + (t40 * t56 * t66 - t21) * t63) + t25 * mrSges(6,2) + t36 * t29;
t14 = m(5) * (-t63 * t21 - t97) + t41 * mrSges(5,2) - t42 * mrSges(5,1) + (t37 * t66 - t38 * t71) * t94 + t76;
t9 = m(5) * (-g(3) * t92 + t81) - t41 * mrSges(5,3) + t50 * mrSges(5,1) - t39 * t84 + t51 * t38 + t65 * t12 + t70 * t11;
t86 = t64 * t14 + t63 * t91 + t9 * t92;
t6 = m(4) * t90 - t54 * mrSges(4,1) - t55 * mrSges(4,2) + t71 * t10 - t66 * t9;
t5 = m(4) * t80 + t55 * mrSges(4,1) - t54 * mrSges(4,2) - t63 * t14 + (t71 * t9 + t91) * t64;
t4 = m(3) * t89 - t58 * mrSges(3,1) - t59 * mrSges(3,2) - t67 * t5 + t72 * t6;
t3 = m(3) * t79 + t59 * mrSges(3,1) - t58 * mrSges(3,2) + t72 * t5 + t67 * t6;
t2 = m(2) * t78 - t75 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t68 * t3 + t73 * t4;
t1 = m(2) * t82 + qJDD(1) * mrSges(2,1) - t75 * mrSges(2,2) + t73 * t3 + t68 * t4;
t7 = [-m(1) * g(1) - t69 * t1 + t74 * t2, t2, t4, t6, t10, t12; -m(1) * g(2) + t74 * t1 + t69 * t2, t1, t3, t5, t9, t11; (-m(1) + t87) * g(3) + t86, t87 * g(3) + t86, t99 * g(3) + t86, -m(4) * g(3) + t86, t14, t76;];
f_new = t7;
