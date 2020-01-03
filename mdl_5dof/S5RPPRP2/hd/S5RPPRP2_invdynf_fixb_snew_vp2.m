% Calculate vector of cutting forces with Newton-Euler
% S5RPPRP2
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2019-12-31 17:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPPRP2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP2_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP2_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP2_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP2_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP2_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP2_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP2_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:49:07
% EndTime: 2019-12-31 17:49:08
% DurationCPUTime: 0.59s
% Computational Cost: add. (5263->121), mult. (11313->151), div. (0->0), fcn. (6783->8), ass. (0->65)
t66 = qJD(1) ^ 2;
t60 = cos(pkin(8));
t55 = t60 ^ 2;
t58 = sin(pkin(8));
t87 = t58 ^ 2 + t55;
t98 = t87 * mrSges(4,3);
t62 = sin(qJ(4));
t94 = cos(qJ(4));
t97 = t58 * t62 - t60 * t94;
t38 = t97 * qJD(1);
t71 = t94 * t58 + t60 * t62;
t39 = t71 * qJD(1);
t26 = t38 * pkin(4) - t39 * qJ(5);
t65 = qJD(4) ^ 2;
t63 = sin(qJ(1));
t64 = cos(qJ(1));
t78 = t63 * g(1) - t64 * g(2);
t44 = qJDD(1) * pkin(1) + t78;
t76 = -t64 * g(1) - t63 * g(2);
t45 = -t66 * pkin(1) + t76;
t59 = sin(pkin(7));
t61 = cos(pkin(7));
t89 = t59 * t44 + t61 * t45;
t23 = -t66 * pkin(2) + qJDD(1) * qJ(3) + t89;
t86 = pkin(6) * qJDD(1);
t57 = -g(3) + qJDD(2);
t83 = qJD(1) * qJD(3);
t88 = t60 * t57 - 0.2e1 * t58 * t83;
t95 = pkin(3) * t66;
t17 = (t60 * t95 - t23 - t86) * t58 + t88;
t80 = t58 * t57 + (t23 + 0.2e1 * t83) * t60;
t18 = -t55 * t95 + t60 * t86 + t80;
t72 = t94 * t17 - t62 * t18;
t96 = m(6) * (-qJDD(4) * pkin(4) - t65 * qJ(5) + t39 * t26 + qJDD(5) - t72);
t92 = -mrSges(5,3) - mrSges(6,2);
t91 = t62 * t17 + t94 * t18;
t27 = t38 * mrSges(6,1) - t39 * mrSges(6,3);
t90 = -t38 * mrSges(5,1) - t39 * mrSges(5,2) - t27;
t85 = t38 * qJD(4);
t84 = t39 * qJD(4);
t32 = t71 * qJDD(1) - t85;
t34 = -qJD(4) * mrSges(5,2) - t38 * mrSges(5,3);
t37 = -t38 * mrSges(6,2) + qJD(4) * mrSges(6,3);
t10 = m(5) * t72 - t96 + t90 * t39 + t92 * t32 + (mrSges(5,1) + mrSges(6,1)) * qJDD(4) + (t34 + t37) * qJD(4);
t74 = -t60 * mrSges(4,1) + t58 * mrSges(4,2);
t73 = qJDD(1) * mrSges(4,3) + t66 * t74;
t31 = t97 * qJDD(1) + t84;
t35 = qJD(4) * mrSges(5,1) - t39 * mrSges(5,3);
t36 = -qJD(4) * mrSges(6,1) + t39 * mrSges(6,2);
t81 = m(6) * (-t65 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) - t38 * t26 + t91) + qJD(4) * t36 + qJDD(4) * mrSges(6,3);
t9 = m(5) * t91 - qJDD(4) * mrSges(5,2) - qJD(4) * t35 + t92 * t31 + t90 * t38 + t81;
t6 = m(4) * t88 + t62 * t9 + t94 * t10 + (-m(4) * t23 - t73) * t58;
t7 = m(4) * t80 - t62 * t10 + t73 * t60 + t94 * t9;
t82 = m(3) * t57 + t58 * t7 + t60 * t6;
t77 = t61 * t44 - t59 * t45;
t75 = qJDD(3) - t77;
t68 = (-pkin(3) * t60 - pkin(2)) * qJDD(1) + (-t87 * pkin(6) - qJ(3)) * t66 + t75;
t70 = t32 * mrSges(6,3) + t39 * t36 - m(6) * (-0.2e1 * qJD(5) * t39 + (-t32 + t85) * qJ(5) + (t31 + t84) * pkin(4) + t68) - t38 * t37 - t31 * mrSges(6,1);
t69 = m(5) * t68 + t31 * mrSges(5,1) + t32 * mrSges(5,2) + t38 * t34 + t39 * t35 - t70;
t67 = m(4) * (-qJDD(1) * pkin(2) - t66 * qJ(3) + t75) + t69;
t8 = -t67 + (-mrSges(3,2) + t98) * t66 + (mrSges(3,1) - t74) * qJDD(1) + m(3) * t77;
t3 = m(3) * t89 - t66 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t58 * t6 + t60 * t7;
t2 = m(2) * t76 - t66 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t61 * t3 - t59 * t8;
t1 = m(2) * t78 + qJDD(1) * mrSges(2,1) - t66 * mrSges(2,2) + t59 * t3 + t61 * t8;
t4 = [-m(1) * g(1) - t63 * t1 + t64 * t2, t2, t3, t7, t9, -t31 * mrSges(6,2) - t38 * t27 + t81; -m(1) * g(2) + t64 * t1 + t63 * t2, t1, t8, t6, t10, -t70; (-m(1) - m(2)) * g(3) + t82, -m(2) * g(3) + t82, t82, t74 * qJDD(1) - t66 * t98 + t67, t69, -qJDD(4) * mrSges(6,1) + t32 * mrSges(6,2) - qJD(4) * t37 + t39 * t27 + t96;];
f_new = t4;
