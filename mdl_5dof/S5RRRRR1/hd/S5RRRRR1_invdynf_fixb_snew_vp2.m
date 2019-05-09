% Calculate vector of cutting forces with Newton-Euler
% S5RRRRR1
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
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
% Datum: 2019-05-04 19:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRRR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR1_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR1_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_invdynf_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR1_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR1_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR1_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:25:47
% EndTime: 2019-05-04 19:25:52
% DurationCPUTime: 1.81s
% Computational Cost: add. (22045->151), mult. (48280->205), div. (0->0), fcn. (36601->10), ass. (0->83)
t77 = sin(qJ(3));
t78 = sin(qJ(2));
t82 = cos(qJ(3));
t83 = cos(qJ(2));
t54 = (t77 * t78 - t82 * t83) * qJD(1);
t55 = (-t77 * t83 - t78 * t82) * qJD(1);
t72 = qJDD(2) + qJDD(3);
t85 = qJD(1) ^ 2;
t79 = sin(qJ(1));
t84 = cos(qJ(1));
t94 = -t84 * g(1) - t79 * g(2);
t61 = -t85 * pkin(1) + t94;
t96 = t83 * g(3) - t78 * t61;
t49 = (t78 * t83 * t85 + qJDD(2)) * pkin(2) + t96;
t103 = t78 * g(3) + t83 * t61;
t50 = (-t83 ^ 2 * t85 - qJD(2) ^ 2) * pkin(2) + t103;
t97 = t82 * t49 - t77 * t50;
t23 = (t54 * t55 + t72) * pkin(3) + t97;
t104 = t77 * t49 + t82 * t50;
t73 = qJD(2) + qJD(3);
t29 = (-t54 ^ 2 - t73 ^ 2) * pkin(3) + t104;
t76 = sin(qJ(4));
t81 = cos(qJ(4));
t105 = t76 * t23 + t81 * t29;
t102 = qJD(1) * t78;
t101 = qJD(1) * t83;
t100 = qJD(1) * qJD(2);
t99 = t79 * g(1) - t84 * g(2);
t98 = t78 * t100;
t95 = qJDD(1) * pkin(1) + t99;
t62 = -t78 * qJDD(1) - t83 * t100;
t63 = -t83 * qJDD(1) + t98;
t37 = t54 * qJD(3) + t82 * t62 + t77 * t63;
t44 = -t54 * mrSges(4,1) + t55 * mrSges(4,2);
t51 = -t73 * mrSges(4,2) + t54 * mrSges(4,3);
t36 = -t55 * qJD(3) - t77 * t62 + t82 * t63;
t43 = t76 * t54 + t81 * t55;
t20 = -t43 * qJD(4) + t81 * t36 - t76 * t37;
t42 = t81 * t54 - t76 * t55;
t21 = t42 * qJD(4) + t76 * t36 + t81 * t37;
t68 = qJD(4) + t73;
t90 = (-t63 - t98) * pkin(2) + t95;
t88 = t90 + (t55 * t73 - t36) * pkin(3);
t13 = (-t42 * t68 - t21) * pkin(6) + (t43 * t68 - t20) * pkin(4) + t88;
t32 = -t42 * pkin(4) - t43 * pkin(6);
t66 = t68 ^ 2;
t67 = qJDD(4) + t72;
t15 = -t66 * pkin(4) + t67 * pkin(6) + t42 * t32 + t105;
t75 = sin(qJ(5));
t80 = cos(qJ(5));
t33 = -t75 * t43 + t80 * t68;
t17 = t33 * qJD(5) + t80 * t21 + t75 * t67;
t19 = qJDD(5) - t20;
t34 = t80 * t43 + t75 * t68;
t24 = -t33 * mrSges(6,1) + t34 * mrSges(6,2);
t40 = qJD(5) - t42;
t26 = -t40 * mrSges(6,2) + t33 * mrSges(6,3);
t11 = m(6) * (t80 * t13 - t75 * t15) - t17 * mrSges(6,3) + t19 * mrSges(6,1) - t34 * t24 + t40 * t26;
t16 = -t34 * qJD(5) - t75 * t21 + t80 * t67;
t27 = t40 * mrSges(6,1) - t34 * mrSges(6,3);
t12 = m(6) * (t75 * t13 + t80 * t15) + t16 * mrSges(6,3) - t19 * mrSges(6,2) + t33 * t24 - t40 * t27;
t31 = -t42 * mrSges(5,1) + t43 * mrSges(5,2);
t39 = t68 * mrSges(5,1) - t43 * mrSges(5,3);
t7 = m(5) * t105 - t67 * mrSges(5,2) + t20 * mrSges(5,3) - t75 * t11 + t80 * t12 + t42 * t31 - t68 * t39;
t38 = -t68 * mrSges(5,2) + t42 * mrSges(5,3);
t92 = t81 * t23 - t76 * t29;
t89 = m(6) * (-t67 * pkin(4) - t66 * pkin(6) + t43 * t32 - t92) - t16 * mrSges(6,1) + t17 * mrSges(6,2) - t33 * t26 + t34 * t27;
t8 = m(5) * t92 + t67 * mrSges(5,1) - t21 * mrSges(5,3) - t43 * t31 + t68 * t38 - t89;
t4 = m(4) * t97 + t72 * mrSges(4,1) - t37 * mrSges(4,3) - t55 * t44 + t73 * t51 + t76 * t7 + t81 * t8;
t52 = t73 * mrSges(4,1) - t55 * mrSges(4,3);
t5 = m(4) * t104 - t72 * mrSges(4,2) + t36 * mrSges(4,3) + t54 * t44 - t73 * t52 + t81 * t7 - t76 * t8;
t60 = (mrSges(3,1) * t83 - mrSges(3,2) * t78) * qJD(1);
t65 = -qJD(2) * mrSges(3,2) - mrSges(3,3) * t101;
t2 = m(3) * t96 + qJDD(2) * mrSges(3,1) - t62 * mrSges(3,3) + qJD(2) * t65 + t60 * t102 + t82 * t4 + t77 * t5;
t64 = qJD(2) * mrSges(3,1) + mrSges(3,3) * t102;
t3 = m(3) * t103 - qJDD(2) * mrSges(3,2) + t63 * mrSges(3,3) - qJD(2) * t64 - t60 * t101 - t77 * t4 + t82 * t5;
t93 = -t83 * t2 - t78 * t3;
t91 = m(5) * t88 - t20 * mrSges(5,1) + t21 * mrSges(5,2) + t80 * t11 + t75 * t12 - t42 * t38 + t43 * t39;
t87 = m(4) * t90 - t36 * mrSges(4,1) + t37 * mrSges(4,2) - t54 * t51 + t55 * t52 + t91;
t86 = m(3) * t95 - t63 * mrSges(3,1) + t62 * mrSges(3,2) + t65 * t101 - t64 * t102 + t87;
t6 = m(2) * t99 + qJDD(1) * mrSges(2,1) - t85 * mrSges(2,2) + t86;
t1 = m(2) * t94 - t85 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t78 * t2 + t83 * t3;
t9 = [-m(1) * g(1) + t84 * t1 - t79 * t6, t1, t3, t5, t7, t12; -m(1) * g(2) + t79 * t1 + t84 * t6, t6, t2, t4, t8, t11; (-m(1) - m(2)) * g(3) + t93, -m(2) * g(3) + t93, t86, t87, t91, t89;];
f_new  = t9;
