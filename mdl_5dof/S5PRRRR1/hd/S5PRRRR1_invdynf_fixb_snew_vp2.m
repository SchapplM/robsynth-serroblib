% Calculate vector of cutting forces with Newton-Euler
% S5PRRRR1
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2019-12-05 17:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRRRR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(2,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR1_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR1_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_invdynf_fixb_snew_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR1_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR1_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR1_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:02:53
% EndTime: 2019-12-05 17:02:55
% DurationCPUTime: 0.62s
% Computational Cost: add. (5844->111), mult. (11284->148), div. (0->0), fcn. (8081->8), ass. (0->62)
t53 = sin(qJ(4));
t54 = sin(qJ(3));
t57 = cos(qJ(4));
t58 = cos(qJ(3));
t36 = (t53 * t54 - t57 * t58) * qJD(2);
t71 = qJD(2) * qJD(3);
t40 = t54 * qJDD(2) + t58 * t71;
t68 = t54 * t71;
t41 = t58 * qJDD(2) - t68;
t73 = qJD(2) * t54;
t44 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t73;
t72 = qJD(2) * t58;
t45 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t72;
t22 = -t36 * qJD(4) + t57 * t40 + t53 * t41;
t37 = (t53 * t58 + t54 * t57) * qJD(2);
t50 = qJD(3) + qJD(4);
t52 = sin(qJ(5));
t56 = cos(qJ(5));
t31 = -t52 * t37 + t56 * t50;
t49 = qJDD(3) + qJDD(4);
t15 = t31 * qJD(5) + t56 * t22 + t52 * t49;
t60 = qJD(2) ^ 2;
t51 = -g(3) + qJDD(1);
t55 = sin(qJ(2));
t59 = cos(qJ(2));
t43 = -t59 * g(1) + t55 * t51;
t67 = t58 * g(2) - t54 * t43;
t29 = (t54 * t58 * t60 + qJDD(3)) * pkin(2) + t67;
t74 = t54 * g(2) + t58 * t43;
t30 = (-t58 ^ 2 * t60 - qJD(3) ^ 2) * pkin(2) + t74;
t17 = t53 * t29 + t57 * t30;
t32 = t56 * t37 + t52 * t50;
t18 = -t31 * mrSges(6,1) + t32 * mrSges(6,2);
t21 = -t37 * qJD(4) - t53 * t40 + t57 * t41;
t20 = qJDD(5) - t21;
t35 = qJD(5) + t36;
t23 = -t35 * mrSges(6,2) + t31 * mrSges(6,3);
t42 = t55 * g(1) + t59 * t51;
t28 = (-t41 + t68) * pkin(2) - t42;
t12 = m(6) * (-t52 * t17 + t56 * t28) - t15 * mrSges(6,3) + t20 * mrSges(6,1) - t32 * t18 + t35 * t23;
t14 = -t32 * qJD(5) - t52 * t22 + t56 * t49;
t24 = t35 * mrSges(6,1) - t32 * mrSges(6,3);
t13 = m(6) * (t56 * t17 + t52 * t28) + t14 * mrSges(6,3) - t20 * mrSges(6,2) + t31 * t18 - t35 * t24;
t33 = -t50 * mrSges(5,2) - t36 * mrSges(5,3);
t34 = t50 * mrSges(5,1) - t37 * mrSges(5,3);
t63 = -m(5) * t28 + t21 * mrSges(5,1) - t22 * mrSges(5,2) - t56 * t12 - t52 * t13 - t36 * t33 - t37 * t34;
t76 = t41 * mrSges(4,1) - t40 * mrSges(4,2) - (t44 * t54 - t45 * t58) * qJD(2) + t63;
t75 = -m(2) - m(3);
t39 = (-mrSges(4,1) * t58 + mrSges(4,2) * t54) * qJD(2);
t26 = t36 * mrSges(5,1) + t37 * mrSges(5,2);
t8 = m(5) * t17 - t49 * mrSges(5,2) + t21 * mrSges(5,3) - t52 * t12 + t56 * t13 - t36 * t26 - t50 * t34;
t16 = -t57 * t29 + t53 * t30;
t62 = t14 * mrSges(6,1) - t15 * mrSges(6,2) + t31 * t23 - t32 * t24;
t9 = t49 * mrSges(5,1) - t22 * mrSges(5,3) - t37 * t26 + t50 * t33 + (-m(5) - m(6)) * t16 + t62;
t4 = m(4) * t67 + qJDD(3) * mrSges(4,1) - t40 * mrSges(4,3) + qJD(3) * t45 - t39 * t73 + t53 * t8 + t57 * t9;
t5 = m(4) * t74 - qJDD(3) * mrSges(4,2) + t41 * mrSges(4,3) - qJD(3) * t44 + t39 * t72 - t53 * t9 + t57 * t8;
t3 = m(3) * t43 - t60 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t54 * t4 + t58 * t5;
t7 = qJDD(2) * mrSges(3,1) - t60 * mrSges(3,2) + (m(3) + m(4)) * t42 + t76;
t70 = m(2) * t51 + t55 * t3 + t59 * t7;
t69 = t59 * t3 - t55 * t7;
t66 = -t58 * t4 - t54 * t5;
t1 = [(-m(1) - m(2)) * g(1) + t69, -m(2) * g(1) + t69, t3, t5, t8, t13; (-m(1) + t75) * g(2) + t66, t75 * g(2) + t66, t7, t4, t9, t12; -m(1) * g(3) + t70, t70, m(3) * g(2) - t66, -m(4) * t42 - t76, -t63, m(6) * t16 - t62;];
f_new = t1;
