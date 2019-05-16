% Calculate vector of cutting forces with Newton-Euler
% S6RRRRRR10V2
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
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
% Datum: 2019-05-08 20:47
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRRR10V2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(6,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10V2_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR10V2_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR10V2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_invdynf_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR10V2_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR10V2_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR10V2_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 20:31:29
% EndTime: 2019-05-08 20:31:52
% DurationCPUTime: 4.22s
% Computational Cost: add. (60464->183), mult. (119187->245), div. (0->0), fcn. (91402->12), ass. (0->99)
t86 = sin(qJ(3));
t87 = sin(qJ(2));
t92 = cos(qJ(3));
t93 = cos(qJ(2));
t68 = (t86 * t87 - t92 * t93) * qJD(1);
t111 = qJD(1) * qJD(2);
t108 = t87 * t111;
t88 = sin(qJ(1));
t94 = cos(qJ(1));
t110 = t88 * g(1) - t94 * g(2);
t71 = -qJDD(1) * pkin(1) - t110;
t75 = t93 * qJDD(1) - t108;
t101 = t71 + (t108 - t75) * pkin(2);
t69 = (t86 * t93 + t87 * t92) * qJD(1);
t74 = t87 * qJDD(1) + t93 * t111;
t49 = -t69 * qJD(3) - t86 * t74 + t92 * t75;
t50 = -t68 * qJD(3) + t92 * t74 + t86 * t75;
t82 = qJD(2) + qJD(3);
t29 = (t68 * t82 - t50) * pkin(5) + (t69 * t82 - t49) * pkin(3) + t101;
t107 = -t94 * g(1) - t88 * g(2);
t95 = qJD(1) ^ 2;
t73 = -t95 * pkin(1) + t107;
t106 = -t93 * g(3) - t87 * t73;
t60 = (t87 * t93 * t95 + qJDD(2)) * pkin(2) + t106;
t109 = -t87 * g(3) + t93 * t73;
t61 = (-t93 ^ 2 * t95 - qJD(2) ^ 2) * pkin(2) + t109;
t114 = t86 * t60 + t92 * t61;
t56 = t68 * pkin(3) - t69 * pkin(5);
t80 = t82 ^ 2;
t81 = qJDD(2) + qJDD(3);
t34 = -t80 * pkin(3) + t81 * pkin(5) - t68 * t56 + t114;
t85 = sin(qJ(4));
t91 = cos(qJ(4));
t24 = t85 * t29 + t91 * t34;
t104 = t92 * t60 - t86 * t61;
t33 = -t81 * pkin(3) - t80 * pkin(5) + t69 * t56 - t104;
t84 = sin(qJ(5));
t90 = cos(qJ(5));
t115 = t90 * t24 + t84 * t33;
t64 = t91 * t69 + t85 * t82;
t37 = -t64 * qJD(4) - t85 * t50 + t91 * t81;
t36 = qJDD(5) - t37;
t67 = qJD(4) + t68;
t45 = -t84 * t64 + t90 * t67;
t46 = t90 * t64 + t84 * t67;
t17 = (-t45 * t46 + t36) * pkin(6) + t115;
t23 = -t91 * t29 + t85 * t34;
t63 = -t85 * t69 + t91 * t82;
t38 = t63 * qJD(4) + t91 * t50 + t85 * t81;
t48 = qJDD(4) - t49;
t27 = t45 * qJD(5) + t90 * t38 + t84 * t48;
t62 = qJD(5) - t63;
t18 = (-t45 * t62 - t27) * pkin(6) + t23;
t83 = sin(qJ(6));
t89 = cos(qJ(6));
t39 = -t83 * t46 + t89 * t62;
t21 = t39 * qJD(6) + t89 * t27 + t83 * t36;
t26 = -t46 * qJD(5) - t84 * t38 + t90 * t48;
t25 = qJDD(6) - t26;
t40 = t89 * t46 + t83 * t62;
t28 = -t39 * mrSges(7,1) + t40 * mrSges(7,2);
t44 = qJD(6) - t45;
t30 = -t44 * mrSges(7,2) + t39 * mrSges(7,3);
t15 = m(7) * (-t83 * t17 + t89 * t18) - t21 * mrSges(7,3) + t25 * mrSges(7,1) - t40 * t28 + t44 * t30;
t20 = -t40 * qJD(6) - t83 * t27 + t89 * t36;
t31 = t44 * mrSges(7,1) - t40 * mrSges(7,3);
t16 = m(7) * (t89 * t17 + t83 * t18) + t20 * mrSges(7,3) - t25 * mrSges(7,2) + t39 * t28 - t44 * t31;
t35 = -t45 * mrSges(6,1) + t46 * mrSges(6,2);
t42 = t62 * mrSges(6,1) - t46 * mrSges(6,3);
t13 = m(6) * t115 - t36 * mrSges(6,2) + t26 * mrSges(6,3) - t83 * t15 + t89 * t16 + t45 * t35 - t62 * t42;
t105 = -t84 * t24 + t90 * t33;
t41 = -t62 * mrSges(6,2) + t45 * mrSges(6,3);
t99 = m(7) * ((-t46 ^ 2 - t62 ^ 2) * pkin(6) - t105) - t20 * mrSges(7,1) + t21 * mrSges(7,2) - t39 * t30 + t40 * t31;
t14 = m(6) * t105 + t36 * mrSges(6,1) - t27 * mrSges(6,3) - t46 * t35 + t62 * t41 - t99;
t43 = -t63 * mrSges(5,1) + t64 * mrSges(5,2);
t52 = t67 * mrSges(5,1) - t64 * mrSges(5,3);
t10 = m(5) * t24 - t48 * mrSges(5,2) + t37 * mrSges(5,3) + t90 * t13 - t84 * t14 + t63 * t43 - t67 * t52;
t51 = -t67 * mrSges(5,2) + t63 * mrSges(5,3);
t98 = t26 * mrSges(6,1) - t27 * mrSges(6,2) - t89 * t15 - t83 * t16 + t45 * t41 - t46 * t42;
t12 = t48 * mrSges(5,1) - t38 * mrSges(5,3) - t64 * t43 + t67 * t51 + (-m(5) - m(6)) * t23 + t98;
t65 = -t82 * mrSges(4,2) - t68 * mrSges(4,3);
t66 = t82 * mrSges(4,1) - t69 * mrSges(4,3);
t100 = -m(4) * t101 + t49 * mrSges(4,1) - t50 * mrSges(4,2) - t85 * t10 - t91 * t12 - t68 * t65 - t69 * t66;
t113 = qJD(1) * t87;
t76 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t113;
t112 = qJD(1) * t93;
t77 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t112;
t117 = m(3) * t71 - t75 * mrSges(3,1) + t74 * mrSges(3,2) + (t87 * t76 - t93 * t77) * qJD(1) - t100;
t54 = t68 * mrSges(4,1) + t69 * mrSges(4,2);
t7 = m(4) * t114 - t81 * mrSges(4,2) + t49 * mrSges(4,3) + t91 * t10 - t85 * t12 - t68 * t54 - t82 * t66;
t72 = (-mrSges(3,1) * t93 + mrSges(3,2) * t87) * qJD(1);
t96 = m(5) * t33 - t37 * mrSges(5,1) + t38 * mrSges(5,2) + t84 * t13 + t90 * t14 - t63 * t51 + t64 * t52;
t8 = m(4) * t104 + t81 * mrSges(4,1) - t50 * mrSges(4,3) - t69 * t54 + t82 * t65 - t96;
t4 = m(3) * t106 + qJDD(2) * mrSges(3,1) - t74 * mrSges(3,3) + qJD(2) * t77 - t72 * t113 + t86 * t7 + t92 * t8;
t5 = m(3) * t109 - qJDD(2) * mrSges(3,2) + t75 * mrSges(3,3) - qJD(2) * t76 + t72 * t112 + t92 * t7 - t86 * t8;
t116 = t93 * t4 + t87 * t5;
t6 = m(2) * t110 + qJDD(1) * mrSges(2,1) - t95 * mrSges(2,2) - t117;
t1 = m(2) * t107 - t95 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t87 * t4 + t93 * t5;
t2 = [-m(1) * g(1) + t94 * t1 - t88 * t6, t1, t5, t7, t10, t13, t16; -m(1) * g(2) + t88 * t1 + t94 * t6, t6, t4, t8, t12, t14, t15; (-m(1) - m(2)) * g(3) + t116, -m(2) * g(3) + t116, t117, -t100, t96, m(6) * t23 - t98, t99;];
f_new  = t2;
