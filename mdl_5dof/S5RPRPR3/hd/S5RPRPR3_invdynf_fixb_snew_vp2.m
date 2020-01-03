% Calculate vector of cutting forces with Newton-Euler
% S5RPRPR3
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2020-01-03 11:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRPR3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR3_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR3_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR3_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR3_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:36:04
% EndTime: 2020-01-03 11:36:07
% DurationCPUTime: 0.72s
% Computational Cost: add. (8362->99), mult. (11760->136), div. (0->0), fcn. (6518->10), ass. (0->65)
t49 = sin(pkin(9));
t51 = cos(pkin(9));
t84 = (t49 ^ 2 + t51 ^ 2) * mrSges(5,3);
t83 = 2 * qJD(4);
t44 = qJDD(1) + qJDD(3);
t82 = t44 * mrSges(5,3);
t47 = qJD(1) + qJD(3);
t81 = t47 * t49;
t80 = t51 * t47;
t48 = -g(1) + qJDD(2);
t79 = t51 * t48;
t55 = sin(qJ(1));
t58 = cos(qJ(1));
t65 = -t58 * g(2) - t55 * g(3);
t34 = qJDD(1) * pkin(1) + t65;
t59 = qJD(1) ^ 2;
t70 = -t55 * g(2) + t58 * g(3);
t35 = -t59 * pkin(1) + t70;
t50 = sin(pkin(8));
t52 = cos(pkin(8));
t67 = t52 * t34 - t50 * t35;
t22 = qJDD(1) * pkin(2) + t67;
t77 = t50 * t34 + t52 * t35;
t23 = -t59 * pkin(2) + t77;
t54 = sin(qJ(3));
t57 = cos(qJ(3));
t78 = t54 * t22 + t57 * t23;
t75 = qJD(5) * t47;
t43 = t47 ^ 2;
t18 = -t43 * pkin(3) + t44 * qJ(4) + t78;
t64 = -t51 * mrSges(5,1) + t49 * mrSges(5,2);
t29 = t64 * t47;
t53 = sin(qJ(5));
t56 = cos(qJ(5));
t27 = (-t44 * t53 - t56 * t75) * t49;
t28 = (t44 * t56 - t53 * t75) * t49;
t66 = -pkin(4) * t51 - pkin(7) * t49;
t30 = t66 * t47;
t61 = m(6) * (-t79 + (t18 + (t83 + t30) * t47) * t49) - t27 * mrSges(6,1) + t28 * mrSges(6,2);
t37 = qJD(5) - t80;
t73 = mrSges(6,3) * t81;
t24 = -t37 * mrSges(6,2) - t53 * t73;
t25 = t37 * mrSges(6,1) - t56 * t73;
t63 = t53 * t24 + t56 * t25;
t10 = m(5) * t79 + (-m(5) * t18 - t82 + (-0.2e1 * m(5) * qJD(4) - t29 - t63) * t47) * t49 - t61;
t71 = t51 * t18 + t49 * t48 + t80 * t83;
t14 = t30 * t80 + t71;
t68 = t57 * t22 - t54 * t23;
t60 = -t43 * qJ(4) + qJDD(4) - t68;
t15 = (-pkin(3) + t66) * t44 + t60;
t36 = -t51 * t44 + qJDD(5);
t72 = (mrSges(6,1) * t53 + mrSges(6,2) * t56) * t81 ^ 2;
t11 = m(6) * (-t53 * t14 + t56 * t15) - t28 * mrSges(6,3) + t36 * mrSges(6,1) - t56 * t72 + t37 * t24;
t12 = m(6) * (t56 * t14 + t53 * t15) + t27 * mrSges(6,3) - t36 * mrSges(6,2) - t53 * t72 - t37 * t25;
t8 = m(5) * t71 + t56 * t12 - t53 * t11 + (t47 * t29 + t82) * t51;
t74 = m(4) * t48 + t51 * t10 + t49 * t8;
t69 = m(3) * t48 + t74;
t62 = m(5) * (-t44 * pkin(3) + t60) + t56 * t11 + t53 * t12;
t6 = m(4) * t68 + (mrSges(4,1) - t64) * t44 + (-mrSges(4,2) + t84) * t43 - t62;
t5 = m(4) * t78 - t43 * mrSges(4,1) - t44 * mrSges(4,2) - t49 * t10 + t51 * t8;
t4 = m(3) * t77 - t59 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t57 * t5 - t54 * t6;
t3 = m(3) * t67 + qJDD(1) * mrSges(3,1) - t59 * mrSges(3,2) + t54 * t5 + t57 * t6;
t2 = m(2) * t65 + qJDD(1) * mrSges(2,1) - t59 * mrSges(2,2) + t52 * t3 + t50 * t4;
t1 = m(2) * t70 - t59 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t50 * t3 + t52 * t4;
t7 = [(-m(1) - m(2)) * g(1) + t69, t1, t4, t5, t8, t12; -m(1) * g(2) + t55 * t1 + t58 * t2, t2, t3, t6, t10, t11; -m(1) * g(3) - t58 * t1 + t55 * t2, -m(2) * g(1) + t69, t69, t74, -t43 * t84 + t64 * t44 + t62, t63 * t81 + t61;];
f_new = t7;
