% Calculate vector of cutting forces with Newton-Euler
% S5RRPRP3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRPRP3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP3_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP3_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP3_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP3_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP3_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:50:58
% EndTime: 2019-12-31 19:50:59
% DurationCPUTime: 0.66s
% Computational Cost: add. (8175->123), mult. (11313->152), div. (0->0), fcn. (6783->8), ass. (0->68)
t55 = qJD(1) + qJD(2);
t51 = t55 ^ 2;
t58 = cos(pkin(8));
t54 = t58 ^ 2;
t57 = sin(pkin(8));
t85 = t57 ^ 2 + t54;
t98 = t85 * mrSges(4,3);
t59 = sin(qJ(4));
t91 = cos(qJ(4));
t97 = t57 * t59 - t58 * t91;
t96 = -m(2) - m(3);
t52 = qJDD(1) + qJDD(2);
t70 = t91 * t57 + t58 * t59;
t38 = t97 * t55;
t82 = t38 * qJD(4);
t31 = t70 * t52 - t82;
t34 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t38;
t37 = -mrSges(6,2) * t38 + qJD(4) * mrSges(6,3);
t39 = t70 * t55;
t61 = sin(qJ(1));
t63 = cos(qJ(1));
t79 = t61 * g(1) - g(2) * t63;
t44 = qJDD(1) * pkin(1) + t79;
t65 = qJD(1) ^ 2;
t75 = -g(1) * t63 - g(2) * t61;
t45 = -pkin(1) * t65 + t75;
t60 = sin(qJ(2));
t62 = cos(qJ(2));
t86 = t60 * t44 + t62 * t45;
t32 = -pkin(2) * t51 + qJ(3) * t52 + t86;
t84 = qJD(3) * t55;
t78 = -g(3) * t58 - 0.2e1 * t57 * t84;
t92 = pkin(7) * t52;
t93 = pkin(3) * t58;
t17 = (t51 * t93 - t32 - t92) * t57 + t78;
t76 = -g(3) * t57 + (t32 + 0.2e1 * t84) * t58;
t18 = -pkin(3) * t51 * t54 + t58 * t92 + t76;
t71 = t91 * t17 - t59 * t18;
t27 = mrSges(6,1) * t38 - mrSges(6,3) * t39;
t87 = -mrSges(5,1) * t38 - mrSges(5,2) * t39 - t27;
t89 = -mrSges(5,3) - mrSges(6,2);
t26 = pkin(4) * t38 - qJ(5) * t39;
t64 = qJD(4) ^ 2;
t94 = m(6) * (-qJDD(4) * pkin(4) - t64 * qJ(5) + t39 * t26 + qJDD(5) - t71);
t10 = m(5) * t71 - t94 + t87 * t39 + t89 * t31 + (mrSges(5,1) + mrSges(6,1)) * qJDD(4) + (t34 + t37) * qJD(4);
t73 = -t58 * mrSges(4,1) + t57 * mrSges(4,2);
t72 = mrSges(4,3) * t52 + t51 * t73;
t83 = qJD(4) * t39;
t30 = t97 * t52 + t83;
t35 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t39;
t36 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t39;
t88 = t59 * t17 + t91 * t18;
t81 = m(6) * (-pkin(4) * t64 + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) - t26 * t38 + t88) + qJD(4) * t36 + qJDD(4) * mrSges(6,3);
t9 = m(5) * t88 - qJDD(4) * mrSges(5,2) - qJD(4) * t35 + t89 * t30 + t87 * t38 + t81;
t6 = m(4) * t78 + t59 * t9 + t91 * t10 + (-m(4) * t32 - t72) * t57;
t7 = m(4) * t76 - t59 * t10 + t72 * t58 + t91 * t9;
t95 = t57 * t7 + t58 * t6;
t77 = t44 * t62 - t60 * t45;
t74 = qJDD(3) - t77;
t67 = (-pkin(2) - t93) * t52 + (-t85 * pkin(7) - qJ(3)) * t51 + t74;
t69 = t31 * mrSges(6,3) + t39 * t36 - m(6) * (-0.2e1 * qJD(5) * t39 + (-t31 + t82) * qJ(5) + (t30 + t83) * pkin(4) + t67) - t38 * t37 - t30 * mrSges(6,1);
t68 = m(5) * t67 + t30 * mrSges(5,1) + t31 * mrSges(5,2) + t38 * t34 + t39 * t35 - t69;
t66 = m(4) * (-pkin(2) * t52 - qJ(3) * t51 + t74) + t68;
t8 = m(3) * t77 - t66 + (-mrSges(3,2) + t98) * t51 + (mrSges(3,1) - t73) * t52;
t3 = m(3) * t86 - t51 * mrSges(3,1) - t52 * mrSges(3,2) - t57 * t6 + t58 * t7;
t2 = m(2) * t75 - t65 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t62 * t3 - t60 * t8;
t1 = m(2) * t79 + qJDD(1) * mrSges(2,1) - t65 * mrSges(2,2) + t60 * t3 + t62 * t8;
t4 = [-m(1) * g(1) - t1 * t61 + t2 * t63, t2, t3, t7, t9, -t30 * mrSges(6,2) - t38 * t27 + t81; -m(1) * g(2) + t1 * t63 + t2 * t61, t1, t8, t6, t10, -t69; (-m(1) + t96) * g(3) + t95, t96 * g(3) + t95, -m(3) * g(3) + t95, -t51 * t98 + t73 * t52 + t66, t68, -qJDD(4) * mrSges(6,1) + t31 * mrSges(6,2) - qJD(4) * t37 + t39 * t27 + t94;];
f_new = t4;
