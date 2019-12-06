% Calculate vector of cutting forces with Newton-Euler
% S5PPRRR3
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
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PPRRR3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR3_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR3_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR3_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR3_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR3_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:16:26
% EndTime: 2019-12-05 15:16:28
% DurationCPUTime: 0.59s
% Computational Cost: add. (5931->99), mult. (10773->139), div. (0->0), fcn. (7077->10), ass. (0->61)
t54 = sin(qJ(4));
t57 = cos(qJ(4));
t71 = qJD(3) * qJD(4);
t69 = t57 * t71;
t37 = t54 * qJDD(3) + t69;
t38 = t57 * qJDD(3) - t54 * t71;
t73 = qJD(3) * t54;
t41 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t73;
t72 = qJD(3) * t57;
t42 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t72;
t59 = qJD(3) ^ 2;
t53 = sin(qJ(5));
t56 = cos(qJ(5));
t34 = (t53 * t57 + t54 * t56) * qJD(3);
t19 = -t34 * qJD(5) - t53 * t37 + t56 * t38;
t33 = (-t53 * t54 + t56 * t57) * qJD(3);
t20 = t33 * qJD(5) + t56 * t37 + t53 * t38;
t46 = qJD(4) + qJD(5);
t29 = -t46 * mrSges(6,2) + t33 * mrSges(6,3);
t30 = t46 * mrSges(6,1) - t34 * mrSges(6,3);
t43 = qJD(4) * pkin(4) - pkin(7) * t73;
t47 = t57 ^ 2;
t50 = sin(pkin(8));
t52 = cos(pkin(8));
t40 = -t52 * g(1) - t50 * g(2);
t48 = -g(3) + qJDD(1);
t49 = sin(pkin(9));
t51 = cos(pkin(9));
t32 = t51 * t40 + t49 * t48;
t66 = t50 * g(1) - t52 * g(2);
t39 = qJDD(2) - t66;
t55 = sin(qJ(3));
t58 = cos(qJ(3));
t67 = -t55 * t32 + t58 * t39;
t63 = -qJDD(3) * pkin(3) - t67;
t61 = t19 * mrSges(6,1) + t33 * t29 - m(6) * (t43 * t73 - t38 * pkin(4) + (-pkin(7) * t47 - pkin(6)) * t59 + t63) - t20 * mrSges(6,2) - t34 * t30;
t76 = (t54 * t41 - t57 * t42) * qJD(3) + m(5) * (-t59 * pkin(6) + t63) - t38 * mrSges(5,1) + t37 * mrSges(5,2) - t61;
t74 = t58 * t32 + t55 * t39;
t22 = -t59 * pkin(3) + qJDD(3) * pkin(6) + t74;
t31 = t49 * t40 - t51 * t48;
t75 = t57 * t22 + t54 * t31;
t10 = m(4) * t67 + qJDD(3) * mrSges(4,1) - t59 * mrSges(4,2) - t76;
t68 = -t54 * t22 + t57 * t31;
t13 = (-t37 + t69) * pkin(7) + (t54 * t57 * t59 + qJDD(4)) * pkin(4) + t68;
t14 = -t47 * t59 * pkin(4) + t38 * pkin(7) - qJD(4) * t43 + t75;
t24 = -t33 * mrSges(6,1) + t34 * mrSges(6,2);
t45 = qJDD(4) + qJDD(5);
t11 = m(6) * (t56 * t13 - t53 * t14) - t20 * mrSges(6,3) + t45 * mrSges(6,1) - t34 * t24 + t46 * t29;
t12 = m(6) * (t53 * t13 + t56 * t14) + t19 * mrSges(6,3) - t45 * mrSges(6,2) + t33 * t24 - t46 * t30;
t36 = (-mrSges(5,1) * t57 + mrSges(5,2) * t54) * qJD(3);
t8 = m(5) * t68 + qJDD(4) * mrSges(5,1) - t37 * mrSges(5,3) + qJD(4) * t42 + t56 * t11 + t53 * t12 - t36 * t73;
t9 = m(5) * t75 - qJDD(4) * mrSges(5,2) + t38 * mrSges(5,3) - qJD(4) * t41 - t53 * t11 + t56 * t12 + t36 * t72;
t6 = m(4) * t74 - t59 * mrSges(4,1) - qJDD(3) * mrSges(4,2) - t54 * t8 + t57 * t9;
t4 = m(3) * t32 - t55 * t10 + t58 * t6;
t65 = -t54 * t9 - t57 * t8;
t7 = (-m(3) - m(4)) * t31 + t65;
t70 = m(2) * t48 + t49 * t4 + t51 * t7;
t62 = m(3) * t39 + t58 * t10 + t55 * t6;
t3 = m(2) * t66 - t62;
t1 = m(2) * t40 + t51 * t4 - t49 * t7;
t2 = [-m(1) * g(1) + t52 * t1 - t50 * t3, t1, t4, t6, t9, t12; -m(1) * g(2) + t50 * t1 + t52 * t3, t3, t7, t10, t8, t11; -m(1) * g(3) + t70, t70, t62, m(4) * t31 - t65, t76, -t61;];
f_new = t2;
