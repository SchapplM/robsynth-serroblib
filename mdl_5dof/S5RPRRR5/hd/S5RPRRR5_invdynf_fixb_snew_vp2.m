% Calculate vector of cutting forces with Newton-Euler
% S5RPRRR5
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
% Datum: 2019-12-05 18:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRRR5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR5_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR5_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR5_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR5_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:16:13
% EndTime: 2019-12-05 18:16:15
% DurationCPUTime: 0.78s
% Computational Cost: add. (11111->112), mult. (15052->152), div. (0->0), fcn. (8503->10), ass. (0->66)
t54 = qJDD(1) + qJDD(3);
t62 = sin(qJ(4));
t66 = cos(qJ(4));
t56 = qJD(1) + qJD(3);
t81 = qJD(4) * t56;
t79 = t66 * t81;
t35 = t62 * t54 + t79;
t36 = t66 * t54 - t62 * t81;
t87 = t56 * t62;
t42 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t87;
t86 = t56 * t66;
t43 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t86;
t52 = t56 ^ 2;
t61 = sin(qJ(5));
t65 = cos(qJ(5));
t33 = (t61 * t66 + t62 * t65) * t56;
t21 = -t33 * qJD(5) - t61 * t35 + t65 * t36;
t32 = (-t61 * t62 + t65 * t66) * t56;
t22 = t32 * qJD(5) + t65 * t35 + t61 * t36;
t55 = qJD(4) + qJD(5);
t30 = -t55 * mrSges(6,2) + t32 * mrSges(6,3);
t31 = t55 * mrSges(6,1) - t33 * mrSges(6,3);
t44 = qJD(4) * pkin(4) - pkin(8) * t87;
t57 = t66 ^ 2;
t64 = sin(qJ(1));
t68 = cos(qJ(1));
t82 = t68 * g(2) + t64 * g(3);
t40 = qJDD(1) * pkin(1) + t82;
t69 = qJD(1) ^ 2;
t78 = t64 * g(2) - t68 * g(3);
t41 = -t69 * pkin(1) + t78;
t59 = sin(pkin(9));
t60 = cos(pkin(9));
t74 = t60 * t40 - t59 * t41;
t28 = qJDD(1) * pkin(2) + t74;
t83 = t59 * t40 + t60 * t41;
t29 = -t69 * pkin(2) + t83;
t63 = sin(qJ(3));
t67 = cos(qJ(3));
t75 = t67 * t28 - t63 * t29;
t72 = -t54 * pkin(3) - t75;
t71 = t21 * mrSges(6,1) + t32 * t30 - m(6) * (t44 * t87 - t36 * pkin(4) + (-pkin(8) * t57 - pkin(7)) * t52 + t72) - t22 * mrSges(6,2) - t33 * t31;
t88 = (t62 * t42 - t66 * t43) * t56 + m(5) * (-t52 * pkin(7) + t72) - t36 * mrSges(5,1) + t35 * mrSges(5,2) - t71;
t84 = t63 * t28 + t67 * t29;
t19 = -t52 * pkin(3) + t54 * pkin(7) + t84;
t58 = -g(1) + qJDD(2);
t85 = t66 * t19 + t62 * t58;
t76 = -t62 * t19 + t66 * t58;
t13 = (-t35 + t79) * pkin(8) + (t52 * t62 * t66 + qJDD(4)) * pkin(4) + t76;
t14 = -t57 * t52 * pkin(4) + t36 * pkin(8) - qJD(4) * t44 + t85;
t24 = -t32 * mrSges(6,1) + t33 * mrSges(6,2);
t53 = qJDD(4) + qJDD(5);
t11 = m(6) * (t65 * t13 - t61 * t14) - t22 * mrSges(6,3) + t53 * mrSges(6,1) - t33 * t24 + t55 * t30;
t12 = m(6) * (t61 * t13 + t65 * t14) + t21 * mrSges(6,3) - t53 * mrSges(6,2) + t32 * t24 - t55 * t31;
t34 = (-mrSges(5,1) * t66 + mrSges(5,2) * t62) * t56;
t8 = m(5) * t76 + qJDD(4) * mrSges(5,1) - t35 * mrSges(5,3) + qJD(4) * t43 + t65 * t11 + t61 * t12 - t34 * t87;
t9 = m(5) * t85 - qJDD(4) * mrSges(5,2) + t36 * mrSges(5,3) - qJD(4) * t42 - t61 * t11 + t65 * t12 + t34 * t86;
t80 = m(4) * t58 + t62 * t9 + t66 * t8;
t77 = m(3) * t58 + t80;
t10 = m(4) * t75 + t54 * mrSges(4,1) - t52 * mrSges(4,2) - t88;
t5 = m(4) * t84 - t52 * mrSges(4,1) - t54 * mrSges(4,2) - t62 * t8 + t66 * t9;
t4 = m(3) * t83 - t69 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t63 * t10 + t67 * t5;
t3 = m(3) * t74 + qJDD(1) * mrSges(3,1) - t69 * mrSges(3,2) + t67 * t10 + t63 * t5;
t2 = m(2) * t82 + qJDD(1) * mrSges(2,1) - t69 * mrSges(2,2) + t60 * t3 + t59 * t4;
t1 = m(2) * t78 - t69 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t59 * t3 + t60 * t4;
t6 = [(-m(1) - m(2)) * g(1) + t77, t1, t4, t5, t9, t12; -m(1) * g(2) - t64 * t1 - t68 * t2, t2, t3, t10, t8, t11; -m(1) * g(3) + t68 * t1 - t64 * t2, -m(2) * g(1) + t77, t77, t80, t88, -t71;];
f_new = t6;
