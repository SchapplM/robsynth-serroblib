% Calculate vector of cutting forces with Newton-Euler
% S5RPPRR7
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPPRR7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR7_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR7_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR7_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR7_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR7_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:59:32
% EndTime: 2019-12-31 17:59:34
% DurationCPUTime: 0.46s
% Computational Cost: add. (3773->110), mult. (6803->140), div. (0->0), fcn. (3520->8), ass. (0->63)
t55 = sin(qJ(1));
t58 = cos(qJ(1));
t73 = t55 * g(1) - g(2) * t58;
t35 = qJDD(1) * pkin(1) + t73;
t60 = qJD(1) ^ 2;
t68 = -g(1) * t58 - g(2) * t55;
t37 = -pkin(1) * t60 + t68;
t51 = sin(pkin(8));
t52 = cos(pkin(8));
t80 = t51 * t35 + t52 * t37;
t85 = -qJDD(1) * qJ(3) - (2 * qJD(3) * qJD(1)) - t80;
t84 = -pkin(2) - pkin(6);
t50 = -g(3) + qJDD(2);
t54 = sin(qJ(4));
t83 = t54 * t50;
t82 = mrSges(3,1) - mrSges(4,2);
t81 = -mrSges(3,2) + mrSges(4,3);
t70 = t52 * t35 - t51 * t37;
t64 = -t60 * qJ(3) + qJDD(3) - t70;
t18 = qJDD(1) * t84 + t64;
t57 = cos(qJ(4));
t79 = t54 * t18 + t57 * t50;
t78 = qJD(1) * t57;
t77 = t54 * qJD(1);
t76 = qJD(1) * qJD(4);
t72 = t54 * t76;
t71 = t57 * t76;
t39 = -qJDD(1) * t54 - t71;
t40 = qJDD(1) * t57 - t72;
t65 = t60 * t84 - t85;
t12 = (-t40 + t72) * pkin(7) + (-t39 + t71) * pkin(4) + t65;
t38 = (pkin(4) * t54 - pkin(7) * t57) * qJD(1);
t59 = qJD(4) ^ 2;
t14 = -pkin(4) * t59 + qJDD(4) * pkin(7) - t38 * t77 + t79;
t53 = sin(qJ(5));
t56 = cos(qJ(5));
t33 = qJD(4) * t56 - t53 * t78;
t22 = qJD(5) * t33 + qJDD(4) * t53 + t40 * t56;
t34 = qJD(4) * t53 + t56 * t78;
t23 = -mrSges(6,1) * t33 + mrSges(6,2) * t34;
t43 = qJD(5) + t77;
t24 = -mrSges(6,2) * t43 + mrSges(6,3) * t33;
t32 = qJDD(5) - t39;
t10 = m(6) * (t12 * t56 - t14 * t53) - t22 * mrSges(6,3) + t32 * mrSges(6,1) - t34 * t23 + t43 * t24;
t21 = -qJD(5) * t34 + qJDD(4) * t56 - t40 * t53;
t25 = mrSges(6,1) * t43 - mrSges(6,3) * t34;
t11 = m(6) * (t12 * t53 + t14 * t56) + t21 * mrSges(6,3) - t32 * mrSges(6,2) + t33 * t23 - t43 * t25;
t36 = (mrSges(5,1) * t54 + mrSges(5,2) * t57) * qJD(1);
t42 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t78;
t6 = m(5) * t79 - qJDD(4) * mrSges(5,2) + t39 * mrSges(5,3) - qJD(4) * t42 - t53 * t10 + t56 * t11 - t36 * t77;
t41 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t77;
t61 = m(6) * (-qJDD(4) * pkin(4) - t59 * pkin(7) + t83 + (qJD(1) * t38 - t18) * t57) - t21 * mrSges(6,1) + t22 * mrSges(6,2) - t33 * t24 + t34 * t25;
t7 = m(5) * (t18 * t57 - t83) - t40 * mrSges(5,3) + qJDD(4) * mrSges(5,1) - t36 * t78 + qJD(4) * t41 - t61;
t69 = m(4) * t50 - t54 * t7 + t57 * t6;
t67 = m(3) * t50 + t69;
t66 = -m(4) * (-qJDD(1) * pkin(2) + t64) - t54 * t6 - t57 * t7;
t63 = m(5) * t65 - t39 * mrSges(5,1) + t40 * mrSges(5,2) + t56 * t10 + t53 * t11 + t41 * t77 + t42 * t78;
t62 = -m(4) * (t60 * pkin(2) + t85) + t63;
t4 = m(3) * t80 + qJDD(1) * t81 - t60 * t82 + t62;
t3 = m(3) * t70 + qJDD(1) * t82 + t60 * t81 + t66;
t2 = m(2) * t68 - t60 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t51 * t3 + t52 * t4;
t1 = m(2) * t73 + qJDD(1) * mrSges(2,1) - t60 * mrSges(2,2) + t52 * t3 + t51 * t4;
t5 = [-m(1) * g(1) - t1 * t55 + t2 * t58, t2, t4, t69, t6, t11; -m(1) * g(2) + t1 * t58 + t2 * t55, t1, t3, -t60 * mrSges(4,2) - qJDD(1) * mrSges(4,3) - t62, t7, t10; (-m(1) - m(2)) * g(3) + t67, -m(2) * g(3) + t67, t67, qJDD(1) * mrSges(4,2) - t60 * mrSges(4,3) - t66, t63, t61;];
f_new = t5;
