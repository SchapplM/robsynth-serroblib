% Calculate vector of cutting forces with Newton-Euler
% S5RRPPR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-05 18:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRPPR2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR2_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR2_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR2_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:20:09
% EndTime: 2019-12-05 18:20:10
% DurationCPUTime: 0.65s
% Computational Cost: add. (8897->100), mult. (11760->136), div. (0->0), fcn. (6518->10), ass. (0->65)
t50 = sin(pkin(9));
t52 = cos(pkin(9));
t85 = (t50 ^ 2 + t52 ^ 2) * mrSges(5,3);
t84 = 2 * qJD(4);
t83 = -m(2) - m(3);
t45 = qJDD(1) + qJDD(2);
t82 = t45 * mrSges(5,3);
t48 = qJD(1) + qJD(2);
t81 = t48 * t50;
t80 = t52 * t48;
t49 = -g(1) + qJDD(3);
t79 = t52 * t49;
t56 = sin(qJ(1));
t59 = cos(qJ(1));
t76 = t59 * g(2) + t56 * g(3);
t34 = qJDD(1) * pkin(1) + t76;
t60 = qJD(1) ^ 2;
t69 = t56 * g(2) - t59 * g(3);
t35 = -t60 * pkin(1) + t69;
t55 = sin(qJ(2));
t58 = cos(qJ(2));
t67 = t58 * t34 - t55 * t35;
t22 = t45 * pkin(2) + t67;
t44 = t48 ^ 2;
t77 = t55 * t34 + t58 * t35;
t23 = -t44 * pkin(2) + t77;
t51 = sin(pkin(8));
t53 = cos(pkin(8));
t78 = t51 * t22 + t53 * t23;
t74 = qJD(5) * t48;
t18 = -t44 * pkin(3) + t45 * qJ(4) + t78;
t65 = -t52 * mrSges(5,1) + t50 * mrSges(5,2);
t29 = t65 * t48;
t54 = sin(qJ(5));
t57 = cos(qJ(5));
t27 = (-t45 * t54 - t57 * t74) * t50;
t28 = (t45 * t57 - t54 * t74) * t50;
t66 = -pkin(4) * t52 - pkin(7) * t50;
t30 = t66 * t48;
t62 = m(6) * (-t79 + (t18 + (t84 + t30) * t48) * t50) - t27 * mrSges(6,1) + t28 * mrSges(6,2);
t37 = qJD(5) - t80;
t72 = mrSges(6,3) * t81;
t24 = -t37 * mrSges(6,2) - t54 * t72;
t25 = t37 * mrSges(6,1) - t57 * t72;
t64 = t54 * t24 + t57 * t25;
t10 = m(5) * t79 + (-m(5) * t18 - t82 + (-0.2e1 * m(5) * qJD(4) - t29 - t64) * t48) * t50 - t62;
t70 = t52 * t18 + t50 * t49 + t80 * t84;
t14 = t30 * t80 + t70;
t68 = t53 * t22 - t51 * t23;
t61 = -t44 * qJ(4) + qJDD(4) - t68;
t15 = (-pkin(3) + t66) * t45 + t61;
t36 = -t52 * t45 + qJDD(5);
t71 = (mrSges(6,1) * t54 + mrSges(6,2) * t57) * t81 ^ 2;
t11 = m(6) * (-t54 * t14 + t57 * t15) - t28 * mrSges(6,3) + t36 * mrSges(6,1) - t57 * t71 + t37 * t24;
t12 = m(6) * (t57 * t14 + t54 * t15) + t27 * mrSges(6,3) - t36 * mrSges(6,2) - t54 * t71 - t37 * t25;
t8 = m(5) * t70 + t57 * t12 - t54 * t11 + (t48 * t29 + t82) * t52;
t73 = m(4) * t49 + t52 * t10 + t50 * t8;
t63 = m(5) * (-t45 * pkin(3) + t61) + t57 * t11 + t54 * t12;
t6 = m(4) * t68 + (mrSges(4,1) - t65) * t45 + (-mrSges(4,2) + t85) * t44 - t63;
t5 = m(4) * t78 - t44 * mrSges(4,1) - t45 * mrSges(4,2) - t50 * t10 + t52 * t8;
t4 = m(3) * t77 - t44 * mrSges(3,1) - t45 * mrSges(3,2) + t53 * t5 - t51 * t6;
t3 = m(3) * t67 + t45 * mrSges(3,1) - t44 * mrSges(3,2) + t51 * t5 + t53 * t6;
t2 = m(2) * t76 + qJDD(1) * mrSges(2,1) - t60 * mrSges(2,2) + t58 * t3 + t55 * t4;
t1 = m(2) * t69 - t60 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t55 * t3 + t58 * t4;
t7 = [(-m(1) + t83) * g(1) + t73, t1, t4, t5, t8, t12; -m(1) * g(2) - t56 * t1 - t59 * t2, t2, t3, t6, t10, t11; -m(1) * g(3) + t59 * t1 - t56 * t2, t83 * g(1) + t73, -m(3) * g(1) + t73, t73, -t44 * t85 + t65 * t45 + t63, t64 * t81 + t62;];
f_new = t7;
