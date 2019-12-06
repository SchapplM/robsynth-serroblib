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
% Datum: 2019-12-05 17:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:51:19
% EndTime: 2019-12-05 17:51:20
% DurationCPUTime: 0.67s
% Computational Cost: add. (8362->99), mult. (11760->136), div. (0->0), fcn. (6518->10), ass. (0->65)
t51 = sin(pkin(9));
t53 = cos(pkin(9));
t86 = (t51 ^ 2 + t53 ^ 2) * mrSges(5,3);
t85 = 2 * qJD(4);
t46 = qJDD(1) + qJDD(3);
t84 = t46 * mrSges(5,3);
t49 = qJD(1) + qJD(3);
t83 = t49 * t51;
t82 = t53 * t49;
t50 = -g(1) + qJDD(2);
t81 = t53 * t50;
t57 = sin(qJ(1));
t60 = cos(qJ(1));
t78 = t60 * g(2) + t57 * g(3);
t34 = qJDD(1) * pkin(1) + t78;
t61 = qJD(1) ^ 2;
t71 = t57 * g(2) - t60 * g(3);
t35 = -t61 * pkin(1) + t71;
t52 = sin(pkin(8));
t54 = cos(pkin(8));
t68 = t54 * t34 - t52 * t35;
t22 = qJDD(1) * pkin(2) + t68;
t79 = t52 * t34 + t54 * t35;
t23 = -t61 * pkin(2) + t79;
t56 = sin(qJ(3));
t59 = cos(qJ(3));
t80 = t56 * t22 + t59 * t23;
t76 = qJD(5) * t49;
t45 = t49 ^ 2;
t18 = -t45 * pkin(3) + t46 * qJ(4) + t80;
t66 = -t53 * mrSges(5,1) + t51 * mrSges(5,2);
t29 = t66 * t49;
t55 = sin(qJ(5));
t58 = cos(qJ(5));
t27 = (-t46 * t55 - t58 * t76) * t51;
t28 = (t46 * t58 - t55 * t76) * t51;
t67 = -pkin(4) * t53 - pkin(7) * t51;
t30 = t67 * t49;
t63 = m(6) * (-t81 + (t18 + (t85 + t30) * t49) * t51) - t27 * mrSges(6,1) + t28 * mrSges(6,2);
t37 = qJD(5) - t82;
t74 = mrSges(6,3) * t83;
t24 = -t37 * mrSges(6,2) - t55 * t74;
t25 = t37 * mrSges(6,1) - t58 * t74;
t65 = t55 * t24 + t58 * t25;
t10 = m(5) * t81 + (-m(5) * t18 - t84 + (-0.2e1 * m(5) * qJD(4) - t29 - t65) * t49) * t51 - t63;
t72 = t53 * t18 + t51 * t50 + t82 * t85;
t14 = t30 * t82 + t72;
t69 = t59 * t22 - t56 * t23;
t62 = -t45 * qJ(4) + qJDD(4) - t69;
t15 = (-pkin(3) + t67) * t46 + t62;
t36 = -t53 * t46 + qJDD(5);
t73 = (mrSges(6,1) * t55 + mrSges(6,2) * t58) * t83 ^ 2;
t11 = m(6) * (-t55 * t14 + t58 * t15) - t28 * mrSges(6,3) + t36 * mrSges(6,1) - t58 * t73 + t37 * t24;
t12 = m(6) * (t58 * t14 + t55 * t15) + t27 * mrSges(6,3) - t36 * mrSges(6,2) - t55 * t73 - t37 * t25;
t8 = m(5) * t72 + t58 * t12 - t55 * t11 + (t49 * t29 + t84) * t53;
t75 = m(4) * t50 + t53 * t10 + t51 * t8;
t70 = m(3) * t50 + t75;
t64 = m(5) * (-t46 * pkin(3) + t62) + t58 * t11 + t55 * t12;
t6 = m(4) * t69 + (mrSges(4,1) - t66) * t46 + (-mrSges(4,2) + t86) * t45 - t64;
t5 = m(4) * t80 - t45 * mrSges(4,1) - t46 * mrSges(4,2) - t51 * t10 + t53 * t8;
t4 = m(3) * t79 - t61 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t59 * t5 - t56 * t6;
t3 = m(3) * t68 + qJDD(1) * mrSges(3,1) - t61 * mrSges(3,2) + t56 * t5 + t59 * t6;
t2 = m(2) * t78 + qJDD(1) * mrSges(2,1) - t61 * mrSges(2,2) + t54 * t3 + t52 * t4;
t1 = m(2) * t71 - t61 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t52 * t3 + t54 * t4;
t7 = [(-m(1) - m(2)) * g(1) + t70, t1, t4, t5, t8, t12; -m(1) * g(2) - t57 * t1 - t60 * t2, t2, t3, t6, t10, t11; -m(1) * g(3) + t60 * t1 - t57 * t2, -m(2) * g(1) + t70, t70, t75, -t45 * t86 + t66 * t46 + t64, t65 * t83 + t63;];
f_new = t7;
