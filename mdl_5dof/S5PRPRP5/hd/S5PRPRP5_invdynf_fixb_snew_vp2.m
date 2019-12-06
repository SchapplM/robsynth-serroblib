% Calculate vector of cutting forces with Newton-Euler
% S5PRPRP5
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRPRP5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP5_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP5_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP5_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP5_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP5_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:37:37
% EndTime: 2019-12-05 15:37:39
% DurationCPUTime: 0.60s
% Computational Cost: add. (4664->115), mult. (10273->145), div. (0->0), fcn. (6585->8), ass. (0->64)
t65 = qJD(2) ^ 2;
t60 = cos(pkin(8));
t55 = t60 ^ 2;
t58 = sin(pkin(8));
t86 = t58 ^ 2 + t55;
t97 = mrSges(4,3) * t86;
t61 = sin(qJ(4));
t93 = cos(qJ(4));
t96 = t58 * t61 - t60 * t93;
t38 = t96 * qJD(2);
t70 = t58 * t93 + t60 * t61;
t39 = t70 * qJD(2);
t23 = t38 * pkin(4) - t39 * qJ(5);
t64 = qJD(4) ^ 2;
t59 = sin(pkin(7));
t85 = cos(pkin(7));
t46 = -g(1) * t85 - t59 * g(2);
t57 = -g(3) + qJDD(1);
t62 = sin(qJ(2));
t63 = cos(qJ(2));
t87 = t63 * t46 + t62 * t57;
t32 = -t65 * pkin(2) + qJDD(2) * qJ(3) + t87;
t84 = pkin(6) * qJDD(2);
t45 = t59 * g(1) - g(2) * t85;
t81 = qJD(2) * qJD(3);
t88 = -t60 * t45 - 0.2e1 * t58 * t81;
t94 = pkin(3) * t65;
t17 = (t60 * t94 - t32 - t84) * t58 + t88;
t78 = -t58 * t45 + (t32 + 0.2e1 * t81) * t60;
t18 = -t55 * t94 + t60 * t84 + t78;
t71 = t17 * t93 - t61 * t18;
t95 = m(6) * (-qJDD(4) * pkin(4) - t64 * qJ(5) + t39 * t23 + qJDD(5) - t71);
t91 = -mrSges(5,3) - mrSges(6,2);
t90 = t61 * t17 + t93 * t18;
t24 = t38 * mrSges(6,1) - t39 * mrSges(6,3);
t89 = -t38 * mrSges(5,1) - t39 * mrSges(5,2) - t24;
t83 = t38 * qJD(4);
t82 = t39 * qJD(4);
t29 = qJDD(2) * t70 - t83;
t34 = -qJD(4) * mrSges(5,2) - t38 * mrSges(5,3);
t37 = -t38 * mrSges(6,2) + qJD(4) * mrSges(6,3);
t10 = m(5) * t71 - t95 + t89 * t39 + t91 * t29 + (mrSges(5,1) + mrSges(6,1)) * qJDD(4) + (t34 + t37) * qJD(4);
t73 = -t60 * mrSges(4,1) + t58 * mrSges(4,2);
t72 = qJDD(2) * mrSges(4,3) + t65 * t73;
t28 = t96 * qJDD(2) + t82;
t35 = qJD(4) * mrSges(5,1) - t39 * mrSges(5,3);
t36 = -qJD(4) * mrSges(6,1) + t39 * mrSges(6,2);
t79 = m(6) * (-t64 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) - t38 * t23 + t90) + qJD(4) * t36 + qJDD(4) * mrSges(6,3);
t9 = m(5) * t90 - qJDD(4) * mrSges(5,2) - qJD(4) * t35 + t28 * t91 + t38 * t89 + t79;
t5 = m(4) * t88 + t61 * t9 + t93 * t10 + (-m(4) * t32 - t72) * t58;
t6 = m(4) * t78 - t61 * t10 + t60 * t72 + t9 * t93;
t3 = m(3) * t87 - t65 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t58 * t5 + t60 * t6;
t76 = -t62 * t46 + t63 * t57;
t74 = qJDD(3) - t76;
t67 = (-pkin(3) * t60 - pkin(2)) * qJDD(2) + (-pkin(6) * t86 - qJ(3)) * t65 + t74;
t69 = t29 * mrSges(6,3) + t39 * t36 - m(6) * (-0.2e1 * qJD(5) * t39 + (-t29 + t83) * qJ(5) + (t28 + t82) * pkin(4) + t67) - t38 * t37 - t28 * mrSges(6,1);
t68 = m(5) * t67 + t28 * mrSges(5,1) + t29 * mrSges(5,2) + t38 * t34 + t39 * t35 - t69;
t66 = m(4) * (-qJDD(2) * pkin(2) - t65 * qJ(3) + t74) + t68;
t8 = -t66 + (-mrSges(3,2) + t97) * t65 + (mrSges(3,1) - t73) * qJDD(2) + m(3) * t76;
t80 = m(2) * t57 + t62 * t3 + t63 * t8;
t75 = -t60 * t5 - t58 * t6;
t4 = (m(2) + m(3)) * t45 + t75;
t1 = m(2) * t46 + t63 * t3 - t62 * t8;
t2 = [-m(1) * g(1) + t1 * t85 - t59 * t4, t1, t3, t6, t9, -t28 * mrSges(6,2) - t38 * t24 + t79; -m(1) * g(2) + t59 * t1 + t4 * t85, t4, t8, t5, t10, -t69; -m(1) * g(3) + t80, t80, -m(3) * t45 - t75, qJDD(2) * t73 - t65 * t97 + t66, t68, -qJDD(4) * mrSges(6,1) + t29 * mrSges(6,2) - qJD(4) * t37 + t39 * t24 + t95;];
f_new = t2;
