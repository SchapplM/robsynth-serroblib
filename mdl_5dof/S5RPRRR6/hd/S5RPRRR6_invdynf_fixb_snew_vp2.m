% Calculate vector of cutting forces with Newton-Euler
% S5RPRRR6
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
% Datum: 2019-12-31 19:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRRR6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR6_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR6_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR6_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR6_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR6_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:01:01
% EndTime: 2019-12-31 19:01:04
% DurationCPUTime: 1.08s
% Computational Cost: add. (12340->138), mult. (24152->186), div. (0->0), fcn. (15329->10), ass. (0->75)
t70 = sin(qJ(4));
t71 = sin(qJ(3));
t74 = cos(qJ(4));
t75 = cos(qJ(3));
t47 = (t70 * t71 - t74 * t75) * qJD(1);
t92 = qJD(1) * qJD(3);
t89 = t75 * t92;
t53 = t71 * qJDD(1) + t89;
t54 = t75 * qJDD(1) - t71 * t92;
t94 = qJD(1) * t71;
t55 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t94;
t93 = qJD(1) * t75;
t56 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t93;
t77 = qJD(1) ^ 2;
t48 = (t70 * t75 + t71 * t74) * qJD(1);
t38 = t47 * pkin(4) - t48 * pkin(8);
t64 = qJD(3) + qJD(4);
t62 = t64 ^ 2;
t63 = qJDD(3) + qJDD(4);
t72 = sin(qJ(1));
t76 = cos(qJ(1));
t90 = t72 * g(1) - t76 * g(2);
t50 = qJDD(1) * pkin(1) + t90;
t86 = -t76 * g(1) - t72 * g(2);
t52 = -t77 * pkin(1) + t86;
t67 = sin(pkin(9));
t68 = cos(pkin(9));
t95 = t67 * t50 + t68 * t52;
t36 = -t77 * pkin(2) + qJDD(1) * pkin(6) + t95;
t66 = -g(3) + qJDD(2);
t88 = -t71 * t36 + t75 * t66;
t22 = (-t53 + t89) * pkin(7) + (t71 * t75 * t77 + qJDD(3)) * pkin(3) + t88;
t57 = qJD(3) * pkin(3) - pkin(7) * t94;
t65 = t75 ^ 2;
t96 = t75 * t36 + t71 * t66;
t23 = -t65 * t77 * pkin(3) + t54 * pkin(7) - qJD(3) * t57 + t96;
t97 = t70 * t22 + t74 * t23;
t16 = -t62 * pkin(4) + t63 * pkin(8) - t47 * t38 + t97;
t29 = -t48 * qJD(4) - t70 * t53 + t74 * t54;
t30 = -t47 * qJD(4) + t74 * t53 + t70 * t54;
t87 = t68 * t50 - t67 * t52;
t82 = -qJDD(1) * pkin(2) - t87;
t79 = -t54 * pkin(3) + t57 * t94 + (-pkin(7) * t65 - pkin(6)) * t77 + t82;
t17 = (t47 * t64 - t30) * pkin(8) + (t48 * t64 - t29) * pkin(4) + t79;
t69 = sin(qJ(5));
t73 = cos(qJ(5));
t39 = -t69 * t48 + t73 * t64;
t19 = t39 * qJD(5) + t73 * t30 + t69 * t63;
t40 = t73 * t48 + t69 * t64;
t26 = -t39 * mrSges(6,1) + t40 * mrSges(6,2);
t28 = qJDD(5) - t29;
t43 = qJD(5) + t47;
t31 = -t43 * mrSges(6,2) + t39 * mrSges(6,3);
t13 = m(6) * (-t69 * t16 + t73 * t17) - t19 * mrSges(6,3) + t28 * mrSges(6,1) - t40 * t26 + t43 * t31;
t18 = -t40 * qJD(5) - t69 * t30 + t73 * t63;
t32 = t43 * mrSges(6,1) - t40 * mrSges(6,3);
t14 = m(6) * (t73 * t16 + t69 * t17) + t18 * mrSges(6,3) - t28 * mrSges(6,2) + t39 * t26 - t43 * t32;
t41 = -t64 * mrSges(5,2) - t47 * mrSges(5,3);
t42 = t64 * mrSges(5,1) - t48 * mrSges(5,3);
t81 = -m(5) * t79 + t29 * mrSges(5,1) - t30 * mrSges(5,2) - t73 * t13 - t69 * t14 - t47 * t41 - t48 * t42;
t98 = (t71 * t55 - t75 * t56) * qJD(1) + m(4) * (-t77 * pkin(6) + t82) - t54 * mrSges(4,1) + t53 * mrSges(4,2) - t81;
t37 = t47 * mrSges(5,1) + t48 * mrSges(5,2);
t85 = t74 * t22 - t70 * t23;
t80 = m(6) * (-t63 * pkin(4) - t62 * pkin(8) + t48 * t38 - t85) - t18 * mrSges(6,1) + t19 * mrSges(6,2) - t39 * t31 + t40 * t32;
t10 = m(5) * t85 + t63 * mrSges(5,1) - t30 * mrSges(5,3) - t48 * t37 + t64 * t41 - t80;
t51 = (-mrSges(4,1) * t75 + mrSges(4,2) * t71) * qJD(1);
t9 = m(5) * t97 - t63 * mrSges(5,2) + t29 * mrSges(5,3) - t69 * t13 + t73 * t14 - t47 * t37 - t64 * t42;
t6 = m(4) * t88 + qJDD(3) * mrSges(4,1) - t53 * mrSges(4,3) + qJD(3) * t56 + t74 * t10 - t51 * t94 + t70 * t9;
t7 = m(4) * t96 - qJDD(3) * mrSges(4,2) + t54 * mrSges(4,3) - qJD(3) * t55 - t70 * t10 + t51 * t93 + t74 * t9;
t91 = m(3) * t66 + t75 * t6 + t71 * t7;
t8 = m(3) * t87 + qJDD(1) * mrSges(3,1) - t77 * mrSges(3,2) - t98;
t3 = m(3) * t95 - t77 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t71 * t6 + t75 * t7;
t2 = m(2) * t86 - t77 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t68 * t3 - t67 * t8;
t1 = m(2) * t90 + qJDD(1) * mrSges(2,1) - t77 * mrSges(2,2) + t67 * t3 + t68 * t8;
t4 = [-m(1) * g(1) - t72 * t1 + t76 * t2, t2, t3, t7, t9, t14; -m(1) * g(2) + t76 * t1 + t72 * t2, t1, t8, t6, t10, t13; (-m(1) - m(2)) * g(3) + t91, -m(2) * g(3) + t91, t91, t98, -t81, t80;];
f_new = t4;
