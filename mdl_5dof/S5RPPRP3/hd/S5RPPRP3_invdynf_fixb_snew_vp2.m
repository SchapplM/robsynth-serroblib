% Calculate vector of cutting forces with Newton-Euler
% S5RPPRP3
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
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
% Datum: 2019-12-31 17:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPPRP3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP3_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP3_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP3_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP3_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP3_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:50:45
% EndTime: 2019-12-31 17:50:46
% DurationCPUTime: 0.35s
% Computational Cost: add. (1909->107), mult. (3431->130), div. (0->0), fcn. (1458->6), ass. (0->52)
t56 = qJD(1) ^ 2;
t52 = sin(qJ(4));
t54 = cos(qJ(4));
t73 = qJD(1) * qJD(4);
t33 = -t52 * qJDD(1) - t54 * t73;
t34 = t54 * qJDD(1) - t52 * t73;
t75 = qJD(1) * t52;
t37 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t75;
t74 = qJD(1) * t54;
t40 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t74;
t53 = sin(qJ(1));
t55 = cos(qJ(1));
t68 = t53 * g(1) - t55 * g(2);
t29 = qJDD(1) * pkin(1) + t68;
t62 = -t55 * g(1) - t53 * g(2);
t32 = -t56 * pkin(1) + t62;
t50 = sin(pkin(7));
t51 = cos(pkin(7));
t77 = t50 * t29 + t51 * t32;
t64 = qJDD(1) * qJ(3) + 0.2e1 * qJD(3) * qJD(1) + t77;
t36 = -qJD(4) * mrSges(6,2) - mrSges(6,3) * t75;
t38 = qJD(4) * pkin(4) - qJ(5) * t74;
t39 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t74;
t48 = t52 ^ 2;
t84 = -pkin(2) - pkin(6);
t65 = m(6) * (t38 * t74 - t33 * pkin(4) + qJDD(5) + (-qJ(5) * t48 + t84) * t56 + t64) + t36 * t75 + t39 * t74 + t34 * mrSges(6,2);
t57 = -(mrSges(5,1) + mrSges(6,1)) * t33 + m(5) * (t84 * t56 + t64) + t37 * t75 + t40 * t74 + t34 * mrSges(5,2) + t65;
t85 = m(4) * (t56 * pkin(2) - t64) - t57;
t82 = pkin(4) * t56;
t69 = -0.2e1 * qJD(1) * qJD(5);
t67 = t51 * t29 - t50 * t32;
t58 = -t56 * qJ(3) + qJDD(3) - t67;
t16 = t84 * qJDD(1) + t58;
t49 = -g(3) + qJDD(2);
t76 = t52 * t16 + t54 * t49;
t81 = t33 * mrSges(6,3) + m(6) * (t33 * qJ(5) - qJD(4) * t38 - t48 * t82 + t52 * t69 + t76);
t80 = mrSges(3,1) - mrSges(4,2);
t78 = -mrSges(3,2) + mrSges(4,3);
t13 = t54 * t16;
t71 = qJD(4) * t36 + qJDD(4) * mrSges(6,1) + m(6) * (t54 * t69 + qJDD(4) * pkin(4) - t34 * qJ(5) + t13 + (-qJ(5) * t73 - t54 * t82 - t49) * t52);
t30 = (mrSges(6,1) * t52 + mrSges(6,2) * t54) * qJD(1);
t66 = qJD(1) * (-t30 - (mrSges(5,1) * t52 + mrSges(5,2) * t54) * qJD(1));
t6 = m(5) * (-t52 * t49 + t13) + qJDD(4) * mrSges(5,1) + qJD(4) * t37 + (-mrSges(5,3) - mrSges(6,3)) * t34 + t54 * t66 + t71;
t7 = m(5) * t76 + t33 * mrSges(5,3) + (-mrSges(5,2) - mrSges(6,2)) * qJDD(4) + (-t40 - t39) * qJD(4) + t52 * t66 + t81;
t63 = m(4) * t49 - t52 * t6 + t54 * t7;
t61 = m(3) * t49 + t63;
t60 = -m(4) * (-qJDD(1) * pkin(2) + t58) - t52 * t7 - t54 * t6;
t4 = m(3) * t77 + t78 * qJDD(1) - t80 * t56 - t85;
t3 = m(3) * t67 + t80 * qJDD(1) + t78 * t56 + t60;
t2 = m(2) * t62 - t56 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t50 * t3 + t51 * t4;
t1 = m(2) * t68 + qJDD(1) * mrSges(2,1) - t56 * mrSges(2,2) + t51 * t3 + t50 * t4;
t5 = [-m(1) * g(1) - t53 * t1 + t55 * t2, t2, t4, t63, t7, -qJDD(4) * mrSges(6,2) - qJD(4) * t39 - t30 * t75 + t81; -m(1) * g(2) + t55 * t1 + t53 * t2, t1, t3, -t56 * mrSges(4,2) - qJDD(1) * mrSges(4,3) + t85, t6, -t34 * mrSges(6,3) - t30 * t74 + t71; (-m(1) - m(2)) * g(3) + t61, -m(2) * g(3) + t61, t61, qJDD(1) * mrSges(4,2) - t56 * mrSges(4,3) - t60, t57, -t33 * mrSges(6,1) + t65;];
f_new = t5;
