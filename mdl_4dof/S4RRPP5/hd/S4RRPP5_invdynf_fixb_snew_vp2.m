% Calculate vector of cutting forces with Newton-Euler
% S4RRPP5
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% Datum: 2019-12-31 17:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4RRPP5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP5_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP5_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP5_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP5_invdynf_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP5_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP5_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP5_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:00:10
% EndTime: 2019-12-31 17:00:11
% DurationCPUTime: 0.37s
% Computational Cost: add. (1057->127), mult. (2172->145), div. (0->0), fcn. (871->4), ass. (0->56)
t54 = qJD(1) ^ 2;
t50 = sin(qJ(1));
t52 = cos(qJ(1));
t66 = t50 * g(1) - t52 * g(2);
t58 = -qJDD(1) * pkin(1) - t66;
t22 = -t54 * pkin(5) + t58;
t49 = sin(qJ(2));
t51 = cos(qJ(2));
t69 = qJD(1) * qJD(2);
t63 = t51 * t69;
t32 = t49 * qJDD(1) + t63;
t64 = t49 * t69;
t33 = t51 * qJDD(1) - t64;
t71 = qJD(1) * t49;
t35 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t71;
t70 = qJD(1) * t51;
t36 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t70;
t40 = mrSges(5,1) * t70 + qJD(2) * mrSges(5,2);
t39 = -mrSges(4,1) * t70 - qJD(2) * mrSges(4,3);
t37 = pkin(3) * t71 - qJD(2) * qJ(4);
t48 = t51 ^ 2;
t85 = -2 * qJD(4);
t86 = -2 * qJD(3);
t89 = pkin(2) * t64 + t71 * t86;
t83 = m(5) * (-t32 * qJ(3) + (-pkin(3) * t48 - pkin(5)) * t54 + (-pkin(2) - qJ(4)) * t33 + (-t37 * t49 + (-qJ(3) * qJD(2) + t85) * t51) * qJD(1) + t58 + t89) - t33 * mrSges(5,3);
t61 = m(4) * (-t33 * pkin(2) + (-t32 - t63) * qJ(3) + t22 + t89) + t39 * t70 - t32 * mrSges(4,3) + t83;
t38 = mrSges(5,1) * t71 - qJD(2) * mrSges(5,3);
t41 = mrSges(4,1) * t71 + qJD(2) * mrSges(4,2);
t73 = -t38 - t41;
t90 = ((t35 + t73) * t49 - (t36 + t40) * t51) * qJD(1) + m(3) * t22 + (mrSges(3,2) - mrSges(5,2)) * t32 - (mrSges(3,1) - mrSges(4,2)) * t33 + t61;
t28 = (-pkin(2) * t51 - qJ(3) * t49) * qJD(1);
t53 = qJD(2) ^ 2;
t60 = -t52 * g(1) - t50 * g(2);
t23 = -t54 * pkin(1) + qJDD(1) * pkin(5) + t60;
t65 = -t49 * g(3) + t51 * t23;
t55 = -t53 * pkin(2) + qJDD(2) * qJ(3) + t28 * t70 + t65;
t31 = (-mrSges(5,2) * t49 - mrSges(5,3) * t51) * qJD(1);
t62 = m(5) * (-t48 * t54 * qJ(4) + t33 * pkin(3) + qJDD(4) + ((2 * qJD(3)) + t37) * qJD(2) + t55) + t31 * t70 + qJD(2) * t38 + qJDD(2) * mrSges(5,2);
t57 = m(4) * (qJD(2) * t86 - t55) - t62;
t29 = (mrSges(4,2) * t51 - mrSges(4,3) * t49) * qJD(1);
t75 = t29 + (-mrSges(3,1) * t51 + mrSges(3,2) * t49) * qJD(1);
t77 = -mrSges(3,3) - mrSges(4,1);
t5 = m(3) * t65 + (-mrSges(3,2) + mrSges(4,3)) * qJDD(2) + (-t35 + t41) * qJD(2) + t75 * t70 + (mrSges(5,1) - t77) * t33 - t57;
t76 = -t51 * g(3) - t49 * t23;
t16 = -qJDD(2) * pkin(2) - t53 * qJ(3) + t28 * t71 + qJDD(3) - t76;
t68 = t31 * t71 + t32 * mrSges(5,1) + m(5) * (qJD(2) * t85 + (-t49 * t51 * t54 - qJDD(2)) * qJ(4) + (t32 - t63) * pkin(3) + t16);
t59 = m(4) * t16 + t68;
t72 = t39 - t40;
t78 = mrSges(4,2) - mrSges(5,3);
t6 = m(3) * t76 + t77 * t32 - t75 * t71 + (mrSges(3,1) - t78) * qJDD(2) + (t36 - t72) * qJD(2) - t59;
t84 = t49 * t5 + t51 * t6;
t82 = t32 * mrSges(5,2);
t81 = t51 * t40;
t2 = m(2) * t66 + qJDD(1) * mrSges(2,1) - t54 * mrSges(2,2) - t90;
t1 = m(2) * t60 - t54 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t49 * t6 + t51 * t5;
t3 = [-m(1) * g(1) + t52 * t1 - t50 * t2, t1, t5, t33 * mrSges(4,2) - t82 + (t73 * t49 - t81) * qJD(1) + t61, -t82 + (-t49 * t38 - t81) * qJD(1) + t83; -m(1) * g(2) + t50 * t1 + t52 * t2, t2, t6, -t29 * t70 - qJDD(2) * mrSges(4,3) - qJD(2) * t41 + (-mrSges(4,1) - mrSges(5,1)) * t33 + t57, -qJDD(2) * mrSges(5,3) - qJD(2) * t40 + t68; (-m(1) - m(2)) * g(3) + t84, -m(2) * g(3) + t84, t90, t32 * mrSges(4,1) + t72 * qJD(2) + t78 * qJDD(2) + t29 * t71 + t59, t33 * mrSges(5,1) + t62;];
f_new = t3;
