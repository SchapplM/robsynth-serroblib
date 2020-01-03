% Calculate vector of cutting forces with Newton-Euler
% S5RPPRR6
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-12-31 17:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPPRR6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR6_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR6_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR6_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR6_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR6_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR6_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR6_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:57:44
% EndTime: 2019-12-31 17:57:46
% DurationCPUTime: 0.94s
% Computational Cost: add. (9823->126), mult. (21471->165), div. (0->0), fcn. (14101->10), ass. (0->74)
t72 = qJD(1) ^ 2;
t63 = cos(pkin(9));
t59 = t63 ^ 2;
t61 = sin(pkin(9));
t92 = t61 ^ 2 + t59;
t97 = t92 * mrSges(4,3);
t66 = sin(qJ(4));
t69 = cos(qJ(4));
t78 = t61 * t66 - t63 * t69;
t44 = t78 * qJD(1);
t79 = t61 * t69 + t63 * t66;
t45 = t79 * qJD(1);
t89 = t45 * qJD(4);
t37 = -t78 * qJDD(1) - t89;
t96 = pkin(3) * t72;
t67 = sin(qJ(1));
t70 = cos(qJ(1));
t85 = t67 * g(1) - g(2) * t70;
t50 = qJDD(1) * pkin(1) + t85;
t83 = -g(1) * t70 - g(2) * t67;
t51 = -pkin(1) * t72 + t83;
t62 = sin(pkin(8));
t64 = cos(pkin(8));
t94 = t62 * t50 + t64 * t51;
t31 = -pkin(2) * t72 + qJDD(1) * qJ(3) + t94;
t91 = pkin(6) * qJDD(1);
t60 = -g(3) + qJDD(2);
t88 = qJD(1) * qJD(3);
t93 = t63 * t60 - 0.2e1 * t61 * t88;
t20 = (t63 * t96 - t31 - t91) * t61 + t93;
t86 = t61 * t60 + (t31 + 0.2e1 * t88) * t63;
t23 = -t59 * t96 + t63 * t91 + t86;
t95 = t66 * t20 + t69 * t23;
t90 = t44 * qJD(4);
t33 = mrSges(5,1) * t44 + mrSges(5,2) * t45;
t38 = t79 * qJDD(1) - t90;
t41 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t44;
t65 = sin(qJ(5));
t68 = cos(qJ(5));
t40 = qJD(4) * t65 + t45 * t68;
t21 = -t40 * qJD(5) + qJDD(4) * t68 - t38 * t65;
t39 = qJD(4) * t68 - t45 * t65;
t22 = t39 * qJD(5) + qJDD(4) * t65 + t38 * t68;
t43 = qJD(5) + t44;
t27 = -mrSges(6,2) * t43 + mrSges(6,3) * t39;
t28 = mrSges(6,1) * t43 - mrSges(6,3) * t40;
t36 = pkin(4) * t44 - pkin(7) * t45;
t71 = qJD(4) ^ 2;
t80 = t20 * t69 - t23 * t66;
t74 = m(6) * (-qJDD(4) * pkin(4) - pkin(7) * t71 + t45 * t36 - t80) - t21 * mrSges(6,1) + t22 * mrSges(6,2) - t39 * t27 + t40 * t28;
t10 = m(5) * t80 + qJDD(4) * mrSges(5,1) - t38 * mrSges(5,3) + qJD(4) * t41 - t45 * t33 - t74;
t81 = -mrSges(4,1) * t63 + mrSges(4,2) * t61;
t77 = mrSges(4,3) * qJDD(1) + t72 * t81;
t16 = -pkin(4) * t71 + qJDD(4) * pkin(7) - t36 * t44 + t95;
t84 = t64 * t50 - t62 * t51;
t82 = qJDD(3) - t84;
t73 = (-pkin(3) * t63 - pkin(2)) * qJDD(1) + (-t92 * pkin(6) - qJ(3)) * t72 + t82;
t17 = (-t38 + t90) * pkin(7) + (-t37 + t89) * pkin(4) + t73;
t26 = -mrSges(6,1) * t39 + mrSges(6,2) * t40;
t35 = qJDD(5) - t37;
t13 = m(6) * (-t16 * t65 + t17 * t68) - t22 * mrSges(6,3) + t35 * mrSges(6,1) - t40 * t26 + t43 * t27;
t14 = m(6) * (t16 * t68 + t17 * t65) + t21 * mrSges(6,3) - t35 * mrSges(6,2) + t39 * t26 - t43 * t28;
t42 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t45;
t9 = m(5) * t95 - qJDD(4) * mrSges(5,2) + t37 * mrSges(5,3) - qJD(4) * t42 - t65 * t13 + t68 * t14 - t44 * t33;
t6 = m(4) * t93 + t66 * t9 + t69 * t10 + (-m(4) * t31 - t77) * t61;
t7 = m(4) * t86 - t66 * t10 + t77 * t63 + t69 * t9;
t87 = m(3) * t60 + t63 * t6 + t61 * t7;
t76 = -m(5) * t73 + t37 * mrSges(5,1) - t38 * mrSges(5,2) - t68 * t13 - t65 * t14 - t44 * t41 - t45 * t42;
t75 = m(4) * (-qJDD(1) * pkin(2) - t72 * qJ(3) + t82) - t76;
t8 = m(3) * t84 + (-mrSges(3,2) + t97) * t72 + (mrSges(3,1) - t81) * qJDD(1) - t75;
t3 = m(3) * t94 - t72 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t61 * t6 + t63 * t7;
t2 = m(2) * t83 - t72 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t64 * t3 - t62 * t8;
t1 = m(2) * t85 + qJDD(1) * mrSges(2,1) - t72 * mrSges(2,2) + t62 * t3 + t64 * t8;
t4 = [-m(1) * g(1) - t1 * t67 + t2 * t70, t2, t3, t7, t9, t14; -m(1) * g(2) + t1 * t70 + t2 * t67, t1, t8, t6, t10, t13; (-m(1) - m(2)) * g(3) + t87, -m(2) * g(3) + t87, t87, t81 * qJDD(1) - t72 * t97 + t75, -t76, t74;];
f_new = t4;
