% Calculate vector of cutting forces with Newton-Euler
% S4RRPP4
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
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4RRPP4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP4_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP4_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_invdynf_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP4_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP4_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP4_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:58:53
% EndTime: 2019-12-31 16:58:54
% DurationCPUTime: 0.37s
% Computational Cost: add. (1057->124), mult. (2186->147), div. (0->0), fcn. (881->4), ass. (0->53)
t53 = sin(qJ(2));
t55 = cos(qJ(2));
t72 = qJD(1) * qJD(2);
t68 = t55 * t72;
t33 = t53 * qJDD(1) + t68;
t34 = t55 * qJDD(1) - t53 * t72;
t74 = qJD(1) * t53;
t39 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t74;
t40 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t74;
t54 = sin(qJ(1));
t56 = cos(qJ(1));
t76 = t54 * g(1) - t56 * g(2);
t64 = qJDD(1) * pkin(1) + t76;
t62 = -t33 * qJ(3) - t64;
t37 = -qJD(2) * pkin(3) - qJ(4) * t74;
t38 = -qJD(2) * mrSges(5,1) - mrSges(5,3) * t74;
t73 = qJD(1) * t55;
t41 = qJD(2) * mrSges(5,2) - mrSges(5,3) * t73;
t52 = t55 ^ 2;
t58 = qJD(1) ^ 2;
t75 = qJ(3) * t55;
t86 = 2 * qJD(3);
t66 = m(5) * (qJDD(4) + (-qJ(4) * t52 + pkin(5)) * t58 + (pkin(2) + pkin(3)) * t34 + (qJD(2) * t75 + (-pkin(2) * qJD(2) + t37 + t86) * t53) * qJD(1) - t62) + t33 * mrSges(5,2) + t34 * mrSges(5,1) + t38 * t74 + t41 * t73;
t84 = t58 * pkin(5);
t63 = m(4) * (-t34 * pkin(2) - t84 + (-0.2e1 * qJD(3) * t53 + (pkin(2) * t53 - t75) * qJD(2)) * qJD(1) + t62) - t34 * mrSges(4,1) - t66;
t43 = mrSges(4,2) * t73 + qJD(2) * mrSges(4,3);
t77 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t73 + t43;
t90 = ((t39 - t40) * t53 - t77 * t55) * qJD(1) + (mrSges(3,2) - mrSges(4,3)) * t33 + m(3) * (-t64 - t84) - t34 * mrSges(3,1) + t63;
t30 = (-mrSges(4,1) * t55 - mrSges(4,3) * t53) * qJD(1);
t29 = (-pkin(2) * t55 - qJ(3) * t53) * qJD(1);
t57 = qJD(2) ^ 2;
t65 = -t56 * g(1) - t54 * g(2);
t22 = -t58 * pkin(1) + qJDD(1) * pkin(5) + t65;
t69 = -t53 * g(3) + t55 * t22;
t60 = -t57 * pkin(2) + qJDD(2) * qJ(3) + qJD(2) * t86 + t29 * t73 + t69;
t88 = m(4) * t60 + qJDD(2) * mrSges(4,3) + qJD(2) * t40 + t30 * t73;
t71 = -0.2e1 * qJD(1) * qJD(4);
t67 = qJD(2) * t38 + qJDD(2) * mrSges(5,2) + m(5) * (-t52 * t58 * pkin(3) - t34 * qJ(4) + qJD(2) * t37 + t55 * t71 + t60) - t34 * mrSges(5,3);
t31 = (mrSges(5,1) * t55 + mrSges(5,2) * t53) * qJD(1);
t79 = -t31 + (-mrSges(3,1) * t55 + mrSges(3,2) * t53) * qJD(1);
t81 = mrSges(3,3) + mrSges(4,2);
t5 = m(3) * t69 - qJDD(2) * mrSges(3,2) - qJD(2) * t39 + t34 * t81 + t73 * t79 + t67 + t88;
t80 = -t55 * g(3) - t53 * t22;
t15 = -qJDD(2) * pkin(2) - t57 * qJ(3) + t29 * t74 + qJDD(3) - t80;
t10 = m(5) * (t53 * t71 + (-t33 + t68) * qJ(4) + (-t53 * t55 * t58 - qJDD(2)) * pkin(3) + t15);
t70 = m(4) * t15 + t10;
t83 = -mrSges(4,1) - mrSges(5,1);
t6 = m(3) * t80 + (mrSges(5,3) - t81) * t33 + (mrSges(3,1) - t83) * qJDD(2) + (t41 + t77) * qJD(2) + (-t30 - t79) * t74 - t70;
t85 = t53 * t5 + t55 * t6;
t61 = -t31 * t73 + t67;
t2 = m(2) * t76 + qJDD(1) * mrSges(2,1) - t58 * mrSges(2,2) - t90;
t1 = m(2) * t65 - t58 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t55 * t5 - t53 * t6;
t3 = [-m(1) * g(1) + t56 * t1 - t54 * t2, t1, t5, t34 * mrSges(4,2) + t61 + t88, t61; -m(1) * g(2) + t54 * t1 + t56 * t2, t2, t6, -t33 * mrSges(4,3) + (-t53 * t40 - t55 * t43) * qJD(1) + t63, -qJDD(2) * mrSges(5,1) - t33 * mrSges(5,3) - qJD(2) * t41 - t31 * t74 + t10; (-m(1) - m(2)) * g(3) + t85, -m(2) * g(3) + t85, t90, (mrSges(4,2) - mrSges(5,3)) * t33 + t83 * qJDD(2) + (-t41 - t43) * qJD(2) + (t30 - t31) * t74 + t70, t66;];
f_new = t3;
