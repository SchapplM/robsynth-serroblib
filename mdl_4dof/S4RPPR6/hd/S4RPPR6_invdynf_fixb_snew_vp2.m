% Calculate vector of cutting forces with Newton-Euler
% S4RPPR6
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% f_new [3x5]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4RPPR6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR6_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR6_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_invdynf_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR6_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR6_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR6_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:38
% EndTime: 2019-12-31 16:40:39
% DurationCPUTime: 0.33s
% Computational Cost: add. (1561->97), mult. (3648->129), div. (0->0), fcn. (2049->6), ass. (0->57)
t83 = mrSges(4,2) + mrSges(3,3);
t47 = sin(pkin(6));
t48 = cos(pkin(6));
t82 = (mrSges(3,2) - mrSges(4,3)) * t47 - (mrSges(3,1) + mrSges(4,1)) * t48;
t63 = -t48 * mrSges(4,1) - t47 * mrSges(4,3);
t37 = t63 * qJD(1);
t38 = (-t48 * mrSges(3,1) + t47 * mrSges(3,2)) * qJD(1);
t74 = qJ(3) * t47;
t62 = -pkin(2) * t48 - t74;
t36 = t62 * qJD(1);
t69 = 0.2e1 * qJD(1) * qJD(2);
t72 = qJD(1) * t47;
t53 = qJD(1) ^ 2;
t50 = sin(qJ(1));
t52 = cos(qJ(1));
t64 = -t52 * g(1) - t50 * g(2);
t35 = -t53 * pkin(1) + qJDD(1) * qJ(2) + t64;
t77 = -t48 * g(3) - t47 * t35;
t14 = t36 * t72 + t47 * t69 + qJDD(3) - t77;
t49 = sin(qJ(4));
t51 = cos(qJ(4));
t80 = pkin(3) * t53;
t10 = (-pkin(5) * qJDD(1) - t48 * t80) * t47 + t14;
t46 = t48 ^ 2;
t66 = -t47 * g(3) + (t35 + t69) * t48;
t71 = qJD(1) * t48;
t59 = t36 * t71 + t66;
t70 = qJDD(1) * t48;
t11 = -pkin(5) * t70 - t46 * t80 + t59;
t60 = -t47 * t49 - t48 * t51;
t33 = t60 * qJD(1);
t61 = t47 * t51 - t48 * t49;
t34 = t61 * qJD(1);
t19 = -t33 * mrSges(5,1) + t34 * mrSges(5,2);
t24 = t33 * qJD(4) + t61 * qJDD(1);
t25 = -qJD(4) * mrSges(5,2) + t33 * mrSges(5,3);
t8 = m(5) * (t51 * t10 - t49 * t11) - t24 * mrSges(5,3) + qJDD(4) * mrSges(5,1) - t34 * t19 + qJD(4) * t25;
t23 = -t34 * qJD(4) + qJDD(1) * t60;
t26 = qJD(4) * mrSges(5,1) - t34 * mrSges(5,3);
t9 = m(5) * (t49 * t10 + t51 * t11) + t23 * mrSges(5,3) - qJDD(4) * mrSges(5,2) + t33 * t19 - qJD(4) * t26;
t56 = m(4) * t14 + t49 * t9 + t51 * t8;
t4 = m(3) * t77 + (-t83 * qJDD(1) + (-0.2e1 * m(3) * qJD(2) - t37 - t38) * qJD(1)) * t47 - t56;
t57 = m(4) * t59 + mrSges(4,2) * t70 + t37 * t71 - t49 * t8 + t51 * t9;
t5 = m(3) * t66 + (qJDD(1) * mrSges(3,3) + qJD(1) * t38) * t48 + t57;
t81 = t48 * t4 + t47 * t5;
t76 = t50 * g(1) - t52 * g(2);
t75 = t47 ^ 2 + t46;
t73 = t53 * qJ(2);
t68 = -qJDD(2) + t76;
t67 = t75 * mrSges(4,2);
t55 = -0.2e1 * qJD(3) * t72 - t68;
t65 = m(5) * ((-t75 * pkin(5) + qJ(2)) * t53 + (t74 + pkin(1) + (pkin(2) + pkin(3)) * t48) * qJDD(1) - t55) + t24 * mrSges(5,2) - t23 * mrSges(5,1) + t34 * t26 - t33 * t25;
t58 = m(4) * (-t73 + (-pkin(1) + t62) * qJDD(1) + t55) - t65;
t54 = m(3) * (-qJDD(1) * pkin(1) - t68 - t73) + t58;
t6 = m(2) * t76 + (t83 * t75 - mrSges(2,2)) * t53 + (mrSges(2,1) - t82) * qJDD(1) - t54;
t1 = m(2) * t64 - t53 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t47 * t4 + t48 * t5;
t2 = [-m(1) * g(1) + t52 * t1 - t50 * t6, t1, t5, t57, t9; -m(1) * g(2) + t50 * t1 + t52 * t6, t6, t4, t63 * qJDD(1) - t53 * t67 + t58, t8; (-m(1) - m(2)) * g(3) + t81, -m(2) * g(3) + t81, (-t75 * mrSges(3,3) - t67) * t53 + t82 * qJDD(1) + t54, (qJDD(1) * mrSges(4,2) + qJD(1) * t37) * t47 + t56, t65;];
f_new = t2;
