% Calculate vector of cutting forces with Newton-Euler
% S5RPPPR4
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
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
% Datum: 2019-12-31 17:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPPPR4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR4_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR4_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR4_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR4_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR4_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR4_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR4_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:45:08
% EndTime: 2019-12-31 17:45:09
% DurationCPUTime: 0.46s
% Computational Cost: add. (3486->98), mult. (6783->126), div. (0->0), fcn. (3768->8), ass. (0->58)
t57 = qJD(1) ^ 2;
t54 = sin(qJ(1));
t56 = cos(qJ(1));
t72 = t54 * g(1) - t56 * g(2);
t34 = qJDD(1) * pkin(1) + t72;
t68 = -t56 * g(1) - t54 * g(2);
t35 = -t57 * pkin(1) + t68;
t50 = sin(pkin(7));
t52 = cos(pkin(7));
t70 = t52 * t34 - t50 * t35;
t61 = -t57 * qJ(3) + qJDD(3) - t70;
t82 = -pkin(2) - qJ(4);
t87 = -(2 * qJD(1) * qJD(4)) + t82 * qJDD(1) + t61;
t49 = sin(pkin(8));
t45 = t49 ^ 2;
t51 = cos(pkin(8));
t79 = t51 ^ 2 + t45;
t71 = t79 * mrSges(5,3);
t80 = t50 * t34 + t52 * t35;
t86 = qJDD(1) * qJ(3) + (2 * qJD(3) * qJD(1)) + t80;
t85 = pkin(4) * t57;
t84 = mrSges(3,1) - mrSges(4,2);
t83 = -mrSges(3,2) + mrSges(4,3);
t81 = t87 * t51;
t78 = qJDD(1) * t49;
t77 = qJDD(1) * t51;
t48 = -g(3) + qJDD(2);
t74 = t51 * t48 + t87 * t49;
t53 = sin(qJ(5));
t55 = cos(qJ(5));
t63 = -qJDD(1) * mrSges(5,3) - t57 * (mrSges(5,1) * t49 + mrSges(5,2) * t51);
t10 = -pkin(6) * t77 + (-t51 * t85 - t48) * t49 + t81;
t11 = -pkin(6) * t78 - t45 * t85 + t74;
t67 = -t49 * t55 - t51 * t53;
t28 = t67 * qJD(1);
t66 = -t49 * t53 + t51 * t55;
t29 = t66 * qJD(1);
t22 = -t28 * mrSges(6,1) + t29 * mrSges(6,2);
t25 = t28 * qJD(5) + t66 * qJDD(1);
t26 = -qJD(5) * mrSges(6,2) + t28 * mrSges(6,3);
t8 = m(6) * (t55 * t10 - t53 * t11) - t25 * mrSges(6,3) + qJDD(5) * mrSges(6,1) - t29 * t22 + qJD(5) * t26;
t24 = -t29 * qJD(5) + t67 * qJDD(1);
t27 = qJD(5) * mrSges(6,1) - t29 * mrSges(6,3);
t9 = m(6) * (t53 * t10 + t55 * t11) + t24 * mrSges(6,3) - qJDD(5) * mrSges(6,2) + t28 * t22 - qJD(5) * t27;
t5 = m(5) * (-t49 * t48 + t81) + t53 * t9 + t55 * t8 + t63 * t51;
t6 = m(5) * t74 + t63 * t49 - t53 * t8 + t55 * t9;
t69 = m(4) * t48 - t49 * t5 + t51 * t6;
t65 = m(3) * t48 + t69;
t64 = qJDD(4) + t86;
t62 = -m(4) * (-qJDD(1) * pkin(2) + t61) - t49 * t6 - t51 * t5;
t60 = t24 * mrSges(6,1) + t28 * t26 - m(6) * (pkin(4) * t78 + (-t79 * pkin(6) + t82) * t57 + t64) - t29 * t27 - t25 * mrSges(6,2);
t59 = m(5) * (t82 * t57 + t64) + mrSges(5,1) * t78 + mrSges(5,2) * t77 - t60;
t58 = m(4) * (t57 * pkin(2) - t86) - t59;
t7 = m(3) * t80 + t83 * qJDD(1) + (-t71 - t84) * t57 - t58;
t3 = m(3) * t70 + t84 * qJDD(1) + t83 * t57 + t62;
t2 = m(2) * t68 - t57 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t50 * t3 + t52 * t7;
t1 = m(2) * t72 + qJDD(1) * mrSges(2,1) - t57 * mrSges(2,2) + t52 * t3 + t50 * t7;
t4 = [-m(1) * g(1) - t54 * t1 + t56 * t2, t2, t7, t69, t6, t9; -m(1) * g(2) + t56 * t1 + t54 * t2, t1, t3, -qJDD(1) * mrSges(4,3) + (-mrSges(4,2) + t71) * t57 + t58, t5, t8; (-m(1) - m(2)) * g(3) + t65, -m(2) * g(3) + t65, t65, qJDD(1) * mrSges(4,2) - t57 * mrSges(4,3) - t62, -t57 * t71 + t59, -t60;];
f_new = t4;
