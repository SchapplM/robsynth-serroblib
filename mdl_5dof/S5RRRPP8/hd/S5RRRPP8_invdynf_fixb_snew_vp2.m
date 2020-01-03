% Calculate vector of cutting forces with Newton-Euler
% S5RRRPP8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 21:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRPP8_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP8_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP8_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP8_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP8_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP8_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP8_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:07:32
% EndTime: 2019-12-31 21:07:34
% DurationCPUTime: 0.68s
% Computational Cost: add. (4947->164), mult. (9625->186), div. (0->0), fcn. (5611->6), ass. (0->75)
t109 = cos(qJ(3));
t71 = sin(qJ(3));
t72 = sin(qJ(2));
t97 = qJD(1) * t72;
t57 = t71 * qJD(2) + t109 * t97;
t74 = cos(qJ(2));
t95 = qJD(1) * qJD(2);
t88 = t74 * t95;
t60 = t72 * qJDD(1) + t88;
t32 = t57 * qJD(3) - t109 * qJDD(2) + t71 * t60;
t56 = -t109 * qJD(2) + t71 * t97;
t33 = -t56 * qJD(3) + t71 * qJDD(2) + t109 * t60;
t96 = t74 * qJD(1);
t65 = -qJD(3) + t96;
t45 = t57 * mrSges(5,1) - t65 * mrSges(5,2);
t108 = t56 * t65;
t112 = -2 * qJD(4);
t59 = (-pkin(2) * t74 - pkin(7) * t72) * qJD(1);
t76 = qJD(2) ^ 2;
t77 = qJD(1) ^ 2;
t73 = sin(qJ(1));
t75 = cos(qJ(1));
t87 = -t75 * g(1) - t73 * g(2);
t52 = -t77 * pkin(1) + qJDD(1) * pkin(6) + t87;
t98 = -t74 * g(3) - t72 * t52;
t23 = -qJDD(2) * pkin(2) - t76 * pkin(7) + t59 * t97 - t98;
t78 = (-t33 - t108) * qJ(4) + t23 + (-t65 * pkin(3) + t112) * t57;
t115 = m(5) * (t32 * pkin(3) + t78) - t33 * mrSges(5,3) - t57 * t45;
t91 = t73 * g(1) - t75 * g(2);
t51 = -qJDD(1) * pkin(1) - t77 * pkin(6) - t91;
t89 = t72 * t95;
t61 = t74 * qJDD(1) - t89;
t62 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t97;
t63 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t96;
t38 = -t56 * mrSges(5,2) - t57 * mrSges(5,3);
t102 = -t56 * mrSges(4,1) - t57 * mrSges(4,2) - t38;
t104 = -mrSges(4,3) - mrSges(5,1);
t105 = mrSges(5,2) - mrSges(6,3);
t39 = t65 * mrSges(4,2) - t56 * mrSges(4,3);
t55 = qJDD(3) - t61;
t36 = t56 * pkin(3) - t57 * qJ(4);
t64 = t65 ^ 2;
t20 = (-t60 - t88) * pkin(7) + (-t61 + t89) * pkin(2) + t51;
t90 = -t72 * g(3) + t74 * t52;
t24 = -t76 * pkin(2) + qJDD(2) * pkin(7) + t59 * t96 + t90;
t86 = t109 * t20 - t71 * t24;
t18 = -t55 * pkin(3) - t64 * qJ(4) + t57 * t36 + qJDD(4) - t86;
t111 = 2 * qJD(5);
t35 = -t57 * mrSges(6,2) + t56 * mrSges(6,3);
t94 = t57 * t35 + t33 * mrSges(6,1) + m(6) * (t65 * t111 + (t56 * t57 - t55) * qJ(5) + (t33 - t108) * pkin(4) + t18);
t84 = m(5) * t18 + t94;
t43 = t56 * mrSges(5,1) + t65 * mrSges(5,3);
t44 = -t56 * mrSges(6,1) - t65 * mrSges(6,2);
t99 = -t43 + t44;
t7 = m(4) * t86 + t102 * t57 + t104 * t33 + (-t39 - t99) * t65 + (mrSges(4,1) - t105) * t55 - t84;
t103 = t109 * t24 + t71 * t20;
t40 = -t65 * mrSges(4,1) - t57 * mrSges(4,3);
t80 = -t64 * pkin(3) + t55 * qJ(4) - t56 * t36 + t103;
t41 = t57 * pkin(4) + t65 * qJ(5);
t42 = t57 * mrSges(6,1) + t65 * mrSges(6,3);
t54 = t56 ^ 2;
t92 = t65 * t42 - m(6) * (-t32 * pkin(4) - t54 * qJ(5) + qJDD(5) + (t112 - t41) * t65 + t80) - t55 * mrSges(6,2);
t83 = m(5) * (0.2e1 * qJD(4) * t65 - t80) + t92;
t8 = m(4) * t103 + (t40 - t45) * t65 + (-mrSges(4,2) + mrSges(5,3)) * t55 + (-t35 + t102) * t56 + (-mrSges(6,1) + t104) * t32 - t83;
t114 = m(3) * t51 - t61 * mrSges(3,1) + t60 * mrSges(3,2) + (t62 * t72 - t63 * t74) * qJD(1) + t109 * t7 + t71 * t8;
t93 = m(6) * (-t54 * pkin(4) + t56 * t111 - t57 * t41 + (pkin(3) + qJ(5)) * t32 + t78) + t32 * mrSges(6,3) + t56 * t44;
t113 = m(4) * t23 + (t40 - t42) * t57 + (t39 - t43) * t56 + (mrSges(4,2) - mrSges(6,2)) * t33 + (mrSges(4,1) - mrSges(5,2)) * t32 + t115 + t93;
t58 = (-mrSges(3,1) * t74 + mrSges(3,2) * t72) * qJD(1);
t4 = m(3) * t90 - qJDD(2) * mrSges(3,2) + t61 * mrSges(3,3) - qJD(2) * t62 + t109 * t8 + t58 * t96 - t71 * t7;
t6 = m(3) * t98 + qJDD(2) * mrSges(3,1) - t60 * mrSges(3,3) + qJD(2) * t63 - t58 * t97 - t113;
t110 = t72 * t4 + t74 * t6;
t82 = -t33 * mrSges(6,2) - t57 * t42 + t93;
t2 = m(2) * t91 + qJDD(1) * mrSges(2,1) - t77 * mrSges(2,2) - t114;
t1 = m(2) * t87 - t77 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t74 * t4 - t72 * t6;
t3 = [-m(1) * g(1) + t75 * t1 - t73 * t2, t1, t4, t8, -t32 * mrSges(5,2) - t56 * t43 + t115 + t82, t82; -m(1) * g(2) + t73 * t1 + t75 * t2, t2, t6, t7, -t55 * mrSges(5,3) + t65 * t45 + (t35 + t38) * t56 + (mrSges(5,1) + mrSges(6,1)) * t32 + t83, -t55 * mrSges(6,3) + t65 * t44 + t94; (-m(1) - m(2)) * g(3) + t110, -m(2) * g(3) + t110, t114, t113, t33 * mrSges(5,1) + t105 * t55 + t57 * t38 + t99 * t65 + t84, -t32 * mrSges(6,1) - t56 * t35 - t92;];
f_new = t3;
