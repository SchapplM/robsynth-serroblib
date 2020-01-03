% Calculate vector of cutting forces with Newton-Euler
% S5RPPRR11
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPPRR11_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR11_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR11_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR11_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR11_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR11_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR11_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:05:33
% EndTime: 2019-12-31 18:05:34
% DurationCPUTime: 0.37s
% Computational Cost: add. (2473->112), mult. (4511->134), div. (0->0), fcn. (2107->6), ass. (0->59)
t60 = qJD(1) ^ 2;
t55 = sin(qJ(1));
t58 = cos(qJ(1));
t83 = t55 * g(1) - t58 * g(2);
t29 = -qJDD(1) * pkin(1) - t60 * qJ(2) + qJDD(2) - t83;
t89 = -qJDD(1) * qJ(3) - (2 * qJD(3) * qJD(1)) + t29;
t70 = -t58 * g(1) - t55 * g(2);
t88 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t70;
t87 = -m(3) - m(4);
t54 = sin(qJ(4));
t86 = t54 * g(3);
t85 = mrSges(3,2) - mrSges(4,3);
t84 = -mrSges(4,2) - mrSges(3,3);
t57 = cos(qJ(4));
t82 = qJD(1) * t57;
t81 = t54 * qJD(1);
t80 = -m(2) + t87;
t79 = qJD(1) * qJD(4);
t76 = mrSges(2,1) - t85;
t72 = t57 * t79;
t39 = -t54 * qJDD(1) - t72;
t73 = t54 * t79;
t40 = t57 * qJDD(1) - t73;
t64 = -t60 * pkin(6) - t89;
t12 = (-t40 + t73) * pkin(7) + (-t39 + t72) * pkin(4) + t64;
t38 = (pkin(4) * t54 - pkin(7) * t57) * qJD(1);
t59 = qJD(4) ^ 2;
t62 = qJDD(3) + (-pkin(1) - qJ(3)) * t60 + t88;
t21 = -qJDD(1) * pkin(6) + t62;
t74 = -t57 * g(3) + t54 * t21;
t14 = -t59 * pkin(4) + qJDD(4) * pkin(7) - t38 * t81 + t74;
t53 = sin(qJ(5));
t56 = cos(qJ(5));
t35 = t56 * qJD(4) - t53 * t82;
t16 = t35 * qJD(5) + t53 * qJDD(4) + t56 * t40;
t36 = t53 * qJD(4) + t56 * t82;
t18 = -t35 * mrSges(6,1) + t36 * mrSges(6,2);
t43 = qJD(5) + t81;
t26 = -t43 * mrSges(6,2) + t35 * mrSges(6,3);
t34 = qJDD(5) - t39;
t10 = m(6) * (t56 * t12 - t53 * t14) - t16 * mrSges(6,3) + t34 * mrSges(6,1) - t36 * t18 + t43 * t26;
t15 = -t36 * qJD(5) + t56 * qJDD(4) - t53 * t40;
t27 = t43 * mrSges(6,1) - t36 * mrSges(6,3);
t11 = m(6) * (t53 * t12 + t56 * t14) + t15 * mrSges(6,3) - t34 * mrSges(6,2) + t35 * t18 - t43 * t27;
t37 = (mrSges(5,1) * t54 + mrSges(5,2) * t57) * qJD(1);
t42 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t82;
t5 = m(5) * t74 - qJDD(4) * mrSges(5,2) + t39 * mrSges(5,3) - qJD(4) * t42 - t53 * t10 + t56 * t11 - t37 * t81;
t41 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t81;
t61 = m(6) * (-qJDD(4) * pkin(4) - t59 * pkin(7) - t86 + (qJD(1) * t38 - t21) * t57) - t15 * mrSges(6,1) + t16 * mrSges(6,2) - t35 * t26 + t36 * t27;
t7 = m(5) * (t57 * t21 + t86) - t40 * mrSges(5,3) + qJDD(4) * mrSges(5,1) - t37 * t82 + qJD(4) * t41 - t61;
t75 = t57 * t5 - t54 * t7;
t71 = m(4) * t62 + qJDD(1) * mrSges(4,2) + t54 * t5 + t57 * t7;
t68 = m(3) * (t60 * pkin(1) - t88) - t71;
t66 = m(5) * t64 - t39 * mrSges(5,1) + t40 * mrSges(5,2) + t56 * t10 + t53 * t11 + t41 * t81 + t42 * t82;
t65 = m(4) * t89 - t66;
t63 = m(3) * t29 + t65;
t2 = m(2) * t83 + (-mrSges(2,2) - t84) * t60 + t76 * qJDD(1) - t63;
t1 = m(2) * t70 + (-mrSges(2,2) + mrSges(3,3)) * qJDD(1) - t76 * t60 - t68;
t3 = [-m(1) * g(1) + t58 * t1 - t55 * t2, t1, t87 * g(3) + t75, -m(4) * g(3) + t75, t5, t11; -m(1) * g(2) + t55 * t1 + t58 * t2, t2, -qJDD(1) * mrSges(3,3) - t85 * t60 + t68, -t60 * mrSges(4,2) - qJDD(1) * mrSges(4,3) + t65, t7, t10; (-m(1) + t80) * g(3) + t75, t80 * g(3) + t75, t85 * qJDD(1) + t84 * t60 + t63, -t60 * mrSges(4,3) + t71, t66, t61;];
f_new = t3;
