% Calculate vector of cutting forces with Newton-Euler
% S5PRRRP3
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRRRP3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP3_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP3_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP3_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP3_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP3_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:43:34
% EndTime: 2019-12-05 16:43:37
% DurationCPUTime: 0.78s
% Computational Cost: add. (6180->127), mult. (12525->166), div. (0->0), fcn. (7843->8), ass. (0->64)
t68 = sin(qJ(3));
t71 = cos(qJ(3));
t88 = qJD(2) * qJD(3);
t83 = t71 * t88;
t49 = t68 * qJDD(2) + t83;
t50 = t71 * qJDD(2) - t68 * t88;
t90 = qJD(2) * t68;
t53 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t90;
t89 = qJD(2) * t71;
t54 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t89;
t73 = qJD(2) ^ 2;
t67 = sin(qJ(4));
t70 = cos(qJ(4));
t43 = (t67 * t71 + t68 * t70) * qJD(2);
t26 = -t43 * qJD(4) - t67 * t49 + t70 * t50;
t42 = (-t67 * t68 + t70 * t71) * qJD(2);
t27 = t42 * qJD(4) + t70 * t49 + t67 * t50;
t62 = qJD(3) + qJD(4);
t36 = -t62 * mrSges(6,2) + t42 * mrSges(6,3);
t37 = -t62 * mrSges(5,2) + t42 * mrSges(5,3);
t40 = t62 * mrSges(5,1) - t43 * mrSges(5,3);
t55 = qJD(3) * pkin(3) - pkin(7) * t90;
t63 = t71 ^ 2;
t65 = sin(pkin(8));
t66 = cos(pkin(8));
t51 = t65 * g(1) - t66 * g(2);
t52 = -t66 * g(1) - t65 * g(2);
t69 = sin(qJ(2));
t72 = cos(qJ(2));
t79 = t72 * t51 - t69 * t52;
t76 = -qJDD(2) * pkin(2) - t79;
t74 = -t50 * pkin(3) + t55 * t90 + (-pkin(7) * t63 - pkin(6)) * t73 + t76;
t38 = t62 * pkin(4) - t43 * qJ(5);
t39 = t62 * mrSges(6,1) - t43 * mrSges(6,3);
t41 = t42 ^ 2;
t84 = m(6) * (-t26 * pkin(4) - t41 * qJ(5) + t43 * t38 + qJDD(5) + t74) + t27 * mrSges(6,2) + t43 * t39;
t75 = m(5) * t74 + t27 * mrSges(5,2) - (mrSges(5,1) + mrSges(6,1)) * t26 + t43 * t40 - (t37 + t36) * t42 + t84;
t100 = t75 + (t68 * t53 - t71 * t54) * qJD(2) - t50 * mrSges(4,1) + t49 * mrSges(4,2) + m(4) * (-t73 * pkin(6) + t76);
t91 = t69 * t51 + t72 * t52;
t34 = -t73 * pkin(2) + qJDD(2) * pkin(6) + t91;
t64 = -g(3) + qJDD(1);
t80 = -t68 * t34 + t71 * t64;
t18 = (-t49 + t83) * pkin(7) + (t68 * t71 * t73 + qJDD(3)) * pkin(3) + t80;
t93 = t71 * t34 + t68 * t64;
t19 = -t63 * t73 * pkin(3) + t50 * pkin(7) - qJD(3) * t55 + t93;
t94 = t67 * t18 + t70 * t19;
t32 = -t42 * mrSges(5,1) + t43 * mrSges(5,2);
t61 = qJDD(3) + qJDD(4);
t31 = -t42 * mrSges(6,1) + t43 * mrSges(6,2);
t85 = m(6) * (-t41 * pkin(4) + t26 * qJ(5) + 0.2e1 * qJD(5) * t42 - t62 * t38 + t94) + t42 * t31 + t26 * mrSges(6,3);
t10 = m(5) * t94 + t26 * mrSges(5,3) + t42 * t32 + (-t40 - t39) * t62 + (-mrSges(5,2) - mrSges(6,2)) * t61 + t85;
t48 = (-mrSges(4,1) * t71 + mrSges(4,2) * t68) * qJD(2);
t81 = t70 * t18 - t67 * t19;
t86 = m(6) * (-0.2e1 * qJD(5) * t43 + (t42 * t62 - t27) * qJ(5) + (t42 * t43 + t61) * pkin(4) + t81) + t62 * t36 + t61 * mrSges(6,1);
t9 = m(5) * t81 + t61 * mrSges(5,1) + t62 * t37 + (-t32 - t31) * t43 + (-mrSges(5,3) - mrSges(6,3)) * t27 + t86;
t6 = m(4) * t80 + qJDD(3) * mrSges(4,1) - t49 * mrSges(4,3) + qJD(3) * t54 + t67 * t10 - t48 * t90 + t70 * t9;
t7 = m(4) * t93 - qJDD(3) * mrSges(4,2) + t50 * mrSges(4,3) - qJD(3) * t53 + t70 * t10 + t48 * t89 - t67 * t9;
t87 = m(3) * t64 + t71 * t6 + t68 * t7;
t82 = m(2) * t64 + t87;
t8 = m(3) * t79 + qJDD(2) * mrSges(3,1) - t73 * mrSges(3,2) - t100;
t3 = m(3) * t91 - t73 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t68 * t6 + t71 * t7;
t2 = m(2) * t52 + t72 * t3 - t69 * t8;
t1 = m(2) * t51 + t69 * t3 + t72 * t8;
t4 = [-m(1) * g(1) - t65 * t1 + t66 * t2, t2, t3, t7, t10, -t61 * mrSges(6,2) - t62 * t39 + t85; -m(1) * g(2) + t66 * t1 + t65 * t2, t1, t8, t6, t9, -t27 * mrSges(6,3) - t43 * t31 + t86; -m(1) * g(3) + t82, t82, t87, t100, t75, -t26 * mrSges(6,1) - t42 * t36 + t84;];
f_new = t4;
