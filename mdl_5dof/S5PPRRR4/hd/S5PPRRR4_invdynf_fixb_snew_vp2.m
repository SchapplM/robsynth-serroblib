% Calculate vector of cutting forces with Newton-Euler
% S5PPRRR4
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PPRRR4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR4_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR4_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR4_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR4_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR4_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:18:55
% EndTime: 2019-12-05 15:18:57
% DurationCPUTime: 1.42s
% Computational Cost: add. (18023->105), mult. (32028->151), div. (0->0), fcn. (24959->14), ass. (0->71)
t51 = sin(pkin(10));
t55 = cos(pkin(10));
t44 = -t55 * g(1) - t51 * g(2);
t50 = sin(pkin(11));
t54 = cos(pkin(11));
t43 = t51 * g(1) - t55 * g(2);
t49 = -g(3) + qJDD(1);
t53 = sin(pkin(5));
t57 = cos(pkin(5));
t69 = t43 * t57 + t49 * t53;
t27 = -t50 * t44 + t69 * t54;
t35 = -t53 * t43 + t57 * t49 + qJDD(2);
t52 = sin(pkin(6));
t56 = cos(pkin(6));
t88 = t27 * t56 + t35 * t52;
t59 = sin(qJ(4));
t62 = cos(qJ(4));
t40 = (-pkin(4) * t62 - pkin(9) * t59) * qJD(3);
t64 = qJD(4) ^ 2;
t77 = t62 * qJD(3);
t65 = qJD(3) ^ 2;
t28 = t54 * t44 + t69 * t50;
t60 = sin(qJ(3));
t63 = cos(qJ(3));
t75 = t63 * t28 + t88 * t60;
t21 = -t65 * pkin(3) + qJDD(3) * pkin(8) + t75;
t23 = -t52 * t27 + t56 * t35;
t79 = t62 * t21 + t59 * t23;
t17 = -t64 * pkin(4) + qJDD(4) * pkin(9) + t40 * t77 + t79;
t86 = -t60 * t28 + t88 * t63;
t20 = -qJDD(3) * pkin(3) - t65 * pkin(8) - t86;
t76 = qJD(3) * qJD(4);
t72 = t62 * t76;
t41 = t59 * qJDD(3) + t72;
t73 = t59 * t76;
t42 = t62 * qJDD(3) - t73;
t18 = (-t41 - t72) * pkin(9) + (-t42 + t73) * pkin(4) + t20;
t58 = sin(qJ(5));
t61 = cos(qJ(5));
t78 = qJD(3) * t59;
t37 = t61 * qJD(4) - t58 * t78;
t30 = t37 * qJD(5) + t58 * qJDD(4) + t61 * t41;
t38 = t58 * qJD(4) + t61 * t78;
t31 = -t37 * mrSges(6,1) + t38 * mrSges(6,2);
t47 = qJD(5) - t77;
t33 = -t47 * mrSges(6,2) + t37 * mrSges(6,3);
t36 = qJDD(5) - t42;
t14 = m(6) * (-t58 * t17 + t61 * t18) - t30 * mrSges(6,3) + t36 * mrSges(6,1) - t38 * t31 + t47 * t33;
t29 = -t38 * qJD(5) + t61 * qJDD(4) - t58 * t41;
t34 = t47 * mrSges(6,1) - t38 * mrSges(6,3);
t15 = m(6) * (t61 * t17 + t58 * t18) + t29 * mrSges(6,3) - t36 * mrSges(6,2) + t37 * t31 - t47 * t34;
t39 = (-mrSges(5,1) * t62 + mrSges(5,2) * t59) * qJD(3);
t45 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t78;
t12 = m(5) * t79 - qJDD(4) * mrSges(5,2) + t42 * mrSges(5,3) - qJD(4) * t45 - t58 * t14 + t61 * t15 + t39 * t77;
t46 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t77;
t80 = t62 * t23;
t66 = m(6) * (-qJDD(4) * pkin(4) - t64 * pkin(9) - t80 + (qJD(3) * t40 + t21) * t59) - t29 * mrSges(6,1) + t30 * mrSges(6,2) - t37 * t33 + t38 * t34;
t13 = m(5) * (-t59 * t21 + t80) - t41 * mrSges(5,3) + qJDD(4) * mrSges(5,1) - t39 * t78 + qJD(4) * t46 - t66;
t10 = m(4) * t23 + t59 * t12 + t62 * t13;
t85 = m(5) * t20 - t42 * mrSges(5,1) + t41 * mrSges(5,2) + t61 * t14 + t58 * t15 + (t59 * t45 - t62 * t46) * qJD(3);
t11 = m(4) * t86 + qJDD(3) * mrSges(4,1) - t65 * mrSges(4,2) - t85;
t9 = m(4) * t75 - t65 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t62 * t12 - t59 * t13;
t71 = t11 * t63 + t60 * t9;
t4 = m(3) * t27 - t52 * t10 + t71 * t56;
t8 = m(3) * t28 - t60 * t11 + t63 * t9;
t87 = t4 * t54 + t50 * t8;
t6 = m(3) * t35 + t56 * t10 + t71 * t52;
t74 = m(2) * t49 + t87 * t53 + t57 * t6;
t2 = m(2) * t44 - t50 * t4 + t54 * t8;
t1 = m(2) * t43 - t53 * t6 + t87 * t57;
t3 = [-m(1) * g(1) - t51 * t1 + t55 * t2, t2, t8, t9, t12, t15; -m(1) * g(2) + t55 * t1 + t51 * t2, t1, t4, t11, t13, t14; -m(1) * g(3) + t74, t74, t6, t10, t85, t66;];
f_new = t3;
