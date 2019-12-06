% Calculate vector of cutting forces with Newton-Euler
% S5PRRRR6
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRRRR6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR6_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR6_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR6_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR6_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR6_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR6_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR6_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:09:38
% EndTime: 2019-12-05 17:09:40
% DurationCPUTime: 0.72s
% Computational Cost: add. (9363->107), mult. (11908->146), div. (0->0), fcn. (7335->10), ass. (0->64)
t54 = qJDD(2) + qJDD(3);
t61 = sin(qJ(4));
t65 = cos(qJ(4));
t56 = qJD(2) + qJD(3);
t79 = qJD(4) * t56;
t76 = t65 * t79;
t37 = t61 * t54 + t76;
t38 = t65 * t54 - t61 * t79;
t85 = t56 * t61;
t42 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t85;
t84 = t56 * t65;
t43 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t84;
t52 = t56 ^ 2;
t60 = sin(qJ(5));
t64 = cos(qJ(5));
t35 = (t60 * t65 + t61 * t64) * t56;
t23 = -t35 * qJD(5) - t60 * t37 + t64 * t38;
t34 = (-t60 * t61 + t64 * t65) * t56;
t24 = t34 * qJD(5) + t64 * t37 + t60 * t38;
t55 = qJD(4) + qJD(5);
t27 = -t55 * mrSges(6,2) + t34 * mrSges(6,3);
t28 = t55 * mrSges(6,1) - t35 * mrSges(6,3);
t44 = qJD(4) * pkin(4) - pkin(8) * t85;
t57 = t65 ^ 2;
t59 = sin(pkin(9));
t80 = cos(pkin(9));
t47 = -t80 * g(1) - t59 * g(2);
t58 = -g(3) + qJDD(1);
t63 = sin(qJ(2));
t67 = cos(qJ(2));
t73 = -t63 * t47 + t67 * t58;
t32 = qJDD(2) * pkin(2) + t73;
t68 = qJD(2) ^ 2;
t81 = t67 * t47 + t63 * t58;
t33 = -t68 * pkin(2) + t81;
t62 = sin(qJ(3));
t66 = cos(qJ(3));
t74 = t66 * t32 - t62 * t33;
t71 = -t54 * pkin(3) - t74;
t70 = t23 * mrSges(6,1) + t34 * t27 - m(6) * (t44 * t85 - t38 * pkin(4) + (-pkin(8) * t57 - pkin(7)) * t52 + t71) - t24 * mrSges(6,2) - t35 * t28;
t86 = (t61 * t42 - t65 * t43) * t56 + m(5) * (-t52 * pkin(7) + t71) - t38 * mrSges(5,1) + t37 * mrSges(5,2) - t70;
t82 = t62 * t32 + t66 * t33;
t21 = -t52 * pkin(3) + t54 * pkin(7) + t82;
t46 = t59 * g(1) - t80 * g(2);
t83 = t65 * t21 - t61 * t46;
t75 = -t61 * t21 - t65 * t46;
t15 = (-t37 + t76) * pkin(8) + (t52 * t61 * t65 + qJDD(4)) * pkin(4) + t75;
t16 = -t57 * t52 * pkin(4) + t38 * pkin(8) - qJD(4) * t44 + t83;
t26 = -t34 * mrSges(6,1) + t35 * mrSges(6,2);
t53 = qJDD(4) + qJDD(5);
t13 = m(6) * (t64 * t15 - t60 * t16) - t24 * mrSges(6,3) + t53 * mrSges(6,1) - t35 * t26 + t55 * t27;
t14 = m(6) * (t60 * t15 + t64 * t16) + t23 * mrSges(6,3) - t53 * mrSges(6,2) + t34 * t26 - t55 * t28;
t36 = (-mrSges(5,1) * t65 + mrSges(5,2) * t61) * t56;
t10 = m(5) * t75 + qJDD(4) * mrSges(5,1) - t37 * mrSges(5,3) + qJD(4) * t43 + t64 * t13 + t60 * t14 - t36 * t85;
t11 = m(5) * t83 - qJDD(4) * mrSges(5,2) + t38 * mrSges(5,3) - qJD(4) * t42 - t60 * t13 + t64 * t14 + t36 * t84;
t78 = m(4) * t46 - t65 * t10 - t61 * t11;
t12 = m(4) * t74 + t54 * mrSges(4,1) - t52 * mrSges(4,2) - t86;
t6 = m(4) * t82 - t52 * mrSges(4,1) - t54 * mrSges(4,2) - t61 * t10 + t65 * t11;
t4 = m(3) * t73 + qJDD(2) * mrSges(3,1) - t68 * mrSges(3,2) + t66 * t12 + t62 * t6;
t5 = m(3) * t81 - t68 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t62 * t12 + t66 * t6;
t77 = m(2) * t58 + t67 * t4 + t63 * t5;
t7 = (m(2) + m(3)) * t46 + t78;
t1 = m(2) * t47 - t63 * t4 + t67 * t5;
t2 = [-m(1) * g(1) + t80 * t1 - t59 * t7, t1, t5, t6, t11, t14; -m(1) * g(2) + t59 * t1 + t80 * t7, t7, t4, t12, t10, t13; -m(1) * g(3) + t77, t77, -m(3) * t46 - t78, -t78, t86, -t70;];
f_new = t2;
