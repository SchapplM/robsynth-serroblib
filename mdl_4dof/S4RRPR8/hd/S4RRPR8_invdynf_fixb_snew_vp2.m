% Calculate vector of cutting forces with Newton-Euler
% S4RRPR8
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
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4RRPR8_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR8_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR8_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR8_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_invdynf_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR8_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR8_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR8_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:07:49
% EndTime: 2019-12-31 17:07:50
% DurationCPUTime: 0.48s
% Computational Cost: add. (2207->125), mult. (4578->163), div. (0->0), fcn. (2311->6), ass. (0->60)
t59 = sin(qJ(2));
t62 = cos(qJ(2));
t77 = qJD(1) * qJD(2);
t75 = t62 * t77;
t40 = t59 * qJDD(1) + t75;
t41 = t62 * qJDD(1) - t59 * t77;
t79 = qJD(1) * t59;
t43 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t79;
t44 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t79;
t60 = sin(qJ(1));
t63 = cos(qJ(1));
t81 = t60 * g(1) - t63 * g(2);
t72 = qJDD(1) * pkin(1) + t81;
t68 = -t40 * qJ(3) - t72;
t58 = sin(qJ(4));
t61 = cos(qJ(4));
t33 = (-t58 * t62 + t59 * t61) * qJD(1);
t20 = -t33 * qJD(4) - t58 * t40 - t61 * t41;
t32 = (-t58 * t59 - t61 * t62) * qJD(1);
t21 = t32 * qJD(4) + t61 * t40 - t58 * t41;
t54 = -qJD(2) + qJD(4);
t25 = -t54 * mrSges(5,2) + t32 * mrSges(5,3);
t26 = t54 * mrSges(5,1) - t33 * mrSges(5,3);
t47 = -qJD(2) * pkin(3) - pkin(6) * t79;
t57 = t62 ^ 2;
t65 = qJD(1) ^ 2;
t80 = qJ(3) * t62;
t89 = 2 * qJD(3);
t74 = m(5) * ((-pkin(6) * t57 + pkin(5)) * t65 + (pkin(2) + pkin(3)) * t41 + (qJD(2) * t80 + (-pkin(2) * qJD(2) + t47 + t89) * t59) * qJD(1) - t68) + t21 * mrSges(5,2) - t20 * mrSges(5,1) + t33 * t26 - t32 * t25;
t87 = t65 * pkin(5);
t71 = m(4) * (-t41 * pkin(2) - t87 + (-0.2e1 * qJD(3) * t59 + (pkin(2) * t59 - t80) * qJD(2)) * qJD(1) + t68) - t41 * mrSges(4,1) - t74;
t78 = qJD(1) * t62;
t46 = mrSges(4,2) * t78 + qJD(2) * mrSges(4,3);
t82 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t78 + t46;
t92 = (-t82 * t62 + (t43 - t44) * t59) * qJD(1) + (mrSges(3,2) - mrSges(4,3)) * t40 + m(3) * (-t72 - t87) - t41 * mrSges(3,1) + t71;
t39 = (-mrSges(3,1) * t62 + mrSges(3,2) * t59) * qJD(1);
t38 = (-mrSges(4,1) * t62 - mrSges(4,3) * t59) * qJD(1);
t37 = (-pkin(2) * t62 - qJ(3) * t59) * qJD(1);
t64 = qJD(2) ^ 2;
t73 = -t63 * g(1) - t60 * g(2);
t35 = -t65 * pkin(1) + qJDD(1) * pkin(5) + t73;
t76 = -t59 * g(3) + t62 * t35;
t67 = -t64 * pkin(2) + qJDD(2) * qJ(3) + qJD(2) * t89 + t37 * t78 + t76;
t12 = -t57 * t65 * pkin(3) - t41 * pkin(6) + qJD(2) * t47 + t67;
t84 = -t62 * g(3) - t59 * t35;
t19 = -qJDD(2) * pkin(2) - t64 * qJ(3) + t37 * t79 + qJDD(3) - t84;
t13 = (-t40 + t75) * pkin(6) + (-t59 * t62 * t65 - qJDD(2)) * pkin(3) + t19;
t24 = -t32 * mrSges(5,1) + t33 * mrSges(5,2);
t53 = -qJDD(2) + qJDD(4);
t8 = m(5) * (-t12 * t58 + t13 * t61) - t21 * mrSges(5,3) + t53 * mrSges(5,1) - t33 * t24 + t54 * t25;
t9 = m(5) * (t12 * t61 + t13 * t58) + t20 * mrSges(5,3) - t53 * mrSges(5,2) + t32 * t24 - t54 * t26;
t70 = m(4) * t67 + qJDD(2) * mrSges(4,3) + qJD(2) * t44 + t38 * t78 - t58 * t8 + t61 * t9;
t85 = mrSges(3,3) + mrSges(4,2);
t4 = m(3) * t76 - qJDD(2) * mrSges(3,2) - qJD(2) * t43 + t39 * t78 + t85 * t41 + t70;
t69 = -m(4) * t19 - t58 * t9 - t61 * t8;
t5 = m(3) * t84 - t85 * t40 + (mrSges(3,1) + mrSges(4,1)) * qJDD(2) + t82 * qJD(2) + (-t38 - t39) * t79 + t69;
t88 = t59 * t4 + t62 * t5;
t6 = m(2) * t81 + qJDD(1) * mrSges(2,1) - t65 * mrSges(2,2) - t92;
t1 = m(2) * t73 - t65 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t62 * t4 - t59 * t5;
t2 = [-m(1) * g(1) + t1 * t63 - t6 * t60, t1, t4, t41 * mrSges(4,2) + t70, t9; -m(1) * g(2) + t1 * t60 + t6 * t63, t6, t5, -t40 * mrSges(4,3) + (-t59 * t44 - t62 * t46) * qJD(1) + t71, t8; (-m(1) - m(2)) * g(3) + t88, -m(2) * g(3) + t88, t92, -qJDD(2) * mrSges(4,1) + t40 * mrSges(4,2) - qJD(2) * t46 + t38 * t79 - t69, t74;];
f_new = t2;
