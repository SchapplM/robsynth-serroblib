% Calculate vector of cutting forces with Newton-Euler
% S5RPPPR3
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3]';
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
% Datum: 2019-12-31 17:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPPPR3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR3_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR3_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR3_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR3_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR3_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR3_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR3_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:43:53
% EndTime: 2019-12-31 17:43:54
% DurationCPUTime: 0.49s
% Computational Cost: add. (3591->111), mult. (7558->146), div. (0->0), fcn. (4237->8), ass. (0->66)
t96 = mrSges(5,2) + mrSges(4,3);
t56 = sin(pkin(8));
t58 = cos(pkin(8));
t95 = (mrSges(4,2) - mrSges(5,3)) * t56 - (mrSges(4,1) + mrSges(5,1)) * t58;
t64 = qJD(1) ^ 2;
t94 = pkin(4) * t64;
t61 = sin(qJ(1));
t63 = cos(qJ(1));
t79 = t61 * g(1) - g(2) * t63;
t45 = qJDD(1) * pkin(1) + t79;
t74 = -g(1) * t63 - g(2) * t61;
t46 = -pkin(1) * t64 + t74;
t57 = sin(pkin(7));
t59 = cos(pkin(7));
t91 = t59 * t45 - t57 * t46;
t90 = t57 * t45 + t59 * t46;
t54 = t58 ^ 2;
t89 = t56 ^ 2 + t54;
t88 = qJ(4) * t56;
t87 = t64 * qJ(3);
t86 = qJD(1) * t56;
t85 = qJD(1) * t58;
t84 = qJDD(1) * t58;
t55 = -g(3) + qJDD(2);
t73 = -t58 * mrSges(5,1) - t56 * mrSges(5,3);
t43 = t73 * qJD(1);
t44 = (-t58 * mrSges(4,1) + t56 * mrSges(4,2)) * qJD(1);
t72 = -pkin(3) * t58 - t88;
t42 = t72 * qJD(1);
t23 = -pkin(2) * t64 + qJDD(1) * qJ(3) + t90;
t77 = -t56 * t23 + t58 * t55;
t81 = 0.2e1 * qJD(1) * qJD(3);
t16 = t42 * t86 + t56 * t81 + qJDD(4) - t77;
t12 = (-pkin(6) * qJDD(1) - t58 * t94) * t56 + t16;
t82 = t56 * t55 + (t23 + t81) * t58;
t76 = t42 * t85 + t82;
t13 = -pkin(6) * t84 - t54 * t94 + t76;
t60 = sin(qJ(5));
t62 = cos(qJ(5));
t70 = -t56 * t60 - t58 * t62;
t36 = t70 * qJD(1);
t71 = t56 * t62 - t58 * t60;
t37 = t71 * qJD(1);
t26 = -mrSges(6,1) * t36 + mrSges(6,2) * t37;
t30 = t36 * qJD(5) + qJDD(1) * t71;
t31 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t36;
t10 = m(6) * (t12 * t62 - t13 * t60) - t30 * mrSges(6,3) + qJDD(5) * mrSges(6,1) - t37 * t26 + qJD(5) * t31;
t29 = -t37 * qJD(5) + qJDD(1) * t70;
t32 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t37;
t11 = m(6) * (t12 * t60 + t13 * t62) + t29 * mrSges(6,3) - qJDD(5) * mrSges(6,2) + t36 * t26 - qJD(5) * t32;
t67 = m(5) * t16 + t62 * t10 + t60 * t11;
t6 = m(4) * t77 + (-t96 * qJDD(1) + (-0.2e1 * m(4) * qJD(3) - t43 - t44) * qJD(1)) * t56 - t67;
t68 = m(5) * t76 + mrSges(5,2) * t84 - t60 * t10 + t62 * t11 + t43 * t85;
t7 = m(4) * t82 + (qJDD(1) * mrSges(4,3) + qJD(1) * t44) * t58 + t68;
t83 = m(3) * t55 + t56 * t7 + t58 * t6;
t80 = -qJDD(3) + t91;
t78 = t89 * mrSges(5,2);
t66 = -0.2e1 * qJD(4) * t86 - t80;
t75 = m(6) * ((-pkin(6) * t89 + qJ(3)) * t64 + (t88 + pkin(2) + (pkin(3) + pkin(4)) * t58) * qJDD(1) - t66) + t30 * mrSges(6,2) - t29 * mrSges(6,1) + t37 * t32 - t36 * t31;
t69 = m(5) * (-t87 + (-pkin(2) + t72) * qJDD(1) + t66) - t75;
t65 = m(4) * (-qJDD(1) * pkin(2) - t80 - t87) + t69;
t8 = m(3) * t91 + (t89 * t96 - mrSges(3,2)) * t64 + (mrSges(3,1) - t95) * qJDD(1) - t65;
t3 = m(3) * t90 - t64 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t56 * t6 + t58 * t7;
t2 = m(2) * t74 - t64 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t59 * t3 - t57 * t8;
t1 = m(2) * t79 + qJDD(1) * mrSges(2,1) - t64 * mrSges(2,2) + t57 * t3 + t59 * t8;
t4 = [-m(1) * g(1) - t1 * t61 + t2 * t63, t2, t3, t7, t68, t11; -m(1) * g(2) + t1 * t63 + t2 * t61, t1, t8, t6, qJDD(1) * t73 - t64 * t78 + t69, t10; (-m(1) - m(2)) * g(3) + t83, -m(2) * g(3) + t83, t83, (-mrSges(4,3) * t89 - t78) * t64 + t95 * qJDD(1) + t65, (qJDD(1) * mrSges(5,2) + qJD(1) * t43) * t56 + t67, t75;];
f_new = t4;
