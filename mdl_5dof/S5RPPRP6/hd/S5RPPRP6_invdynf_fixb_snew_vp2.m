% Calculate vector of cutting forces with Newton-Euler
% S5RPPRP6
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
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
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPPRP6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP6_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP6_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP6_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP6_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP6_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:54:59
% EndTime: 2019-12-31 17:55:00
% DurationCPUTime: 0.51s
% Computational Cost: add. (3335->120), mult. (7260->144), div. (0->0), fcn. (4231->6), ass. (0->63)
t61 = qJD(1) ^ 2;
t57 = sin(qJ(1));
t59 = cos(qJ(1));
t78 = t57 * g(1) - t59 * g(2);
t68 = -t61 * qJ(2) + qJDD(2) - t78;
t90 = -pkin(1) - qJ(3);
t99 = -(2 * qJD(1) * qJD(3)) + qJDD(1) * t90 + t68;
t54 = sin(pkin(7));
t50 = t54 ^ 2;
t55 = cos(pkin(7));
t87 = t55 ^ 2 + t50;
t77 = t87 * mrSges(4,3);
t75 = -t59 * g(1) - t57 * g(2);
t98 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t75;
t97 = -m(2) - m(3);
t56 = sin(qJ(4));
t58 = cos(qJ(4));
t73 = t54 * t58 + t55 * t56;
t39 = t73 * qJD(1);
t72 = -t54 * t56 + t55 * t58;
t40 = t72 * qJD(1);
t20 = pkin(4) * t39 - qJ(5) * t40;
t60 = qJD(4) ^ 2;
t81 = t54 * g(3) + t55 * t99;
t95 = pkin(3) * t61;
t14 = (-pkin(6) * qJDD(1) - t54 * t95) * t55 + t81;
t76 = -t55 * g(3) + t54 * t99;
t84 = qJDD(1) * t54;
t15 = -pkin(6) * t84 - t50 * t95 + t76;
t74 = t14 * t58 - t15 * t56;
t96 = m(6) * (-qJDD(4) * pkin(4) - qJ(5) * t60 + t20 * t40 + qJDD(5) - t74);
t94 = mrSges(4,2) * t55;
t93 = mrSges(2,1) - mrSges(3,2);
t92 = -mrSges(2,2) + mrSges(3,3);
t91 = -mrSges(5,3) - mrSges(6,2);
t89 = t14 * t56 + t15 * t58;
t21 = mrSges(6,1) * t39 - mrSges(6,3) * t40;
t88 = -mrSges(5,1) * t39 - mrSges(5,2) * t40 - t21;
t86 = t39 * qJD(4);
t85 = t40 * qJD(4);
t34 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t40;
t82 = m(6) * (-pkin(4) * t60 + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) - t20 * t39 + t89) + qJD(4) * t34 + qJDD(4) * mrSges(6,3);
t25 = qJDD(1) * t73 + t85;
t33 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t40;
t6 = m(5) * t89 - qJDD(4) * mrSges(5,2) - qJD(4) * t33 + t25 * t91 + t39 * t88 + t82;
t26 = qJDD(1) * t72 - t86;
t32 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t39;
t35 = -mrSges(6,2) * t39 + qJD(4) * mrSges(6,3);
t7 = m(5) * t74 - t96 + t88 * t40 + t91 * t26 + (mrSges(5,1) + mrSges(6,1)) * qJDD(4) + (t32 + t35) * qJD(4);
t70 = -qJDD(1) * mrSges(4,3) - t61 * (mrSges(4,1) * t54 + t94);
t3 = m(4) * t81 + t55 * t70 + t56 * t6 + t58 * t7;
t4 = m(4) * t76 + t54 * t70 - t56 * t7 + t58 * t6;
t79 = -t54 * t3 + t4 * t55;
t69 = -m(3) * (-qJDD(1) * pkin(1) + t68) - t55 * t3 - t54 * t4;
t66 = qJDD(3) + t98;
t64 = pkin(3) * t84 + (-pkin(6) * t87 + t90) * t61 + t66;
t67 = t26 * mrSges(6,3) + t40 * t34 - t39 * t35 - t25 * mrSges(6,1) - m(6) * (-0.2e1 * qJD(5) * t40 + (-t26 + t86) * qJ(5) + (t25 + t85) * pkin(4) + t64);
t65 = m(5) * t64 + t25 * mrSges(5,1) + t26 * mrSges(5,2) + t39 * t32 + t40 * t33 - t67;
t63 = m(4) * (t61 * t90 + t66) + mrSges(4,1) * t84 + qJDD(1) * t94 + t65;
t62 = m(3) * (t61 * pkin(1) - t98) - t63;
t5 = -t62 + (-t77 - t93) * t61 + t92 * qJDD(1) + m(2) * t75;
t1 = m(2) * t78 + qJDD(1) * t93 + t61 * t92 + t69;
t2 = [-m(1) * g(1) - t1 * t57 + t5 * t59, t5, -m(3) * g(3) + t79, t4, t6, -t25 * mrSges(6,2) - t39 * t21 + t82; -m(1) * g(2) + t1 * t59 + t5 * t57, t1, -qJDD(1) * mrSges(3,3) + t62 + (-mrSges(3,2) + t77) * t61, t3, t7, -t67; (-m(1) + t97) * g(3) + t79, g(3) * t97 + t79, qJDD(1) * mrSges(3,2) - t61 * mrSges(3,3) - t69, -t61 * t77 + t63, t65, -qJDD(4) * mrSges(6,1) + t26 * mrSges(6,2) - qJD(4) * t35 + t40 * t21 + t96;];
f_new = t2;
