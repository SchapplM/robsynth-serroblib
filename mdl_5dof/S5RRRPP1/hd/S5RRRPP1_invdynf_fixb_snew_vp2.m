% Calculate vector of cutting forces with Newton-Euler
% S5RRRPP1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Datum: 2019-12-31 20:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRPP1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP1_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP1_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP1_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP1_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP1_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:49:27
% EndTime: 2019-12-31 20:49:28
% DurationCPUTime: 0.74s
% Computational Cost: add. (9335->136), mult. (12482->175), div. (0->0), fcn. (7135->8), ass. (0->68)
t57 = qJDD(1) + qJDD(2);
t62 = sin(qJ(3));
t65 = cos(qJ(3));
t58 = qJD(1) + qJD(2);
t83 = qJD(3) * t58;
t44 = t57 * t62 + t65 * t83;
t45 = t57 * t65 - t62 * t83;
t90 = t58 * t62;
t52 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t90;
t89 = t58 * t65;
t53 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t89;
t56 = t58 ^ 2;
t61 = sin(pkin(8));
t84 = cos(pkin(8));
t31 = t44 * t61 - t45 * t84;
t32 = t44 * t84 + t45 * t61;
t39 = t61 * t90 - t84 * t89;
t34 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t39;
t40 = (t61 * t65 + t62 * t84) * t58;
t35 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t40;
t51 = qJD(3) * pkin(3) - qJ(4) * t90;
t60 = t65 ^ 2;
t64 = sin(qJ(1));
t67 = cos(qJ(1));
t80 = t64 * g(1) - g(2) * t67;
t49 = qJDD(1) * pkin(1) + t80;
t69 = qJD(1) ^ 2;
t77 = -g(1) * t67 - g(2) * t64;
t50 = -pkin(1) * t69 + t77;
t63 = sin(qJ(2));
t66 = cos(qJ(2));
t78 = t66 * t49 - t63 * t50;
t75 = -t57 * pkin(2) - t78;
t71 = -t45 * pkin(3) + qJDD(4) + t51 * t90 + (-qJ(4) * t60 - pkin(7)) * t56 + t75;
t36 = -qJD(3) * mrSges(6,1) + mrSges(6,2) * t40;
t37 = -mrSges(6,2) * t39 + qJD(3) * mrSges(6,3);
t73 = t32 * mrSges(6,3) + t40 * t36 - m(6) * (-0.2e1 * qJD(5) * t40 + (qJD(3) * t39 - t32) * qJ(5) + (qJD(3) * t40 + t31) * pkin(4) + t71) - t39 * t37 - t31 * mrSges(6,1);
t72 = m(5) * t71 + t31 * mrSges(5,1) + t32 * mrSges(5,2) + t39 * t34 + t40 * t35 - t73;
t96 = (t52 * t62 - t53 * t65) * t58 + m(4) * (-t56 * pkin(7) + t75) - t45 * mrSges(4,1) + t44 * mrSges(4,2) + t72;
t95 = -2 * qJD(4);
t94 = -m(2) - m(3);
t85 = t63 * t49 + t66 * t50;
t28 = -pkin(2) * t56 + pkin(7) * t57 + t85;
t88 = t62 * t28;
t91 = pkin(3) * t56;
t17 = qJDD(3) * pkin(3) - t44 * qJ(4) - t88 + (qJ(4) * t83 + t62 * t91 - g(3)) * t65;
t79 = -g(3) * t62 + t65 * t28;
t18 = qJ(4) * t45 - qJD(3) * t51 - t60 * t91 + t79;
t74 = t17 * t84 - t61 * t18;
t25 = mrSges(6,1) * t39 - mrSges(6,3) * t40;
t86 = -mrSges(5,1) * t39 - mrSges(5,2) * t40 - t25;
t87 = -mrSges(5,3) - mrSges(6,2);
t24 = pkin(4) * t39 - qJ(5) * t40;
t68 = qJD(3) ^ 2;
t92 = m(6) * (-qJDD(3) * pkin(4) - t68 * qJ(5) + qJDD(5) + ((2 * qJD(4)) + t24) * t40 - t74);
t10 = m(5) * t74 - t92 + (m(5) * t95 + t86) * t40 + t87 * t32 + (mrSges(5,1) + mrSges(6,1)) * qJDD(3) + (t34 + t37) * qJD(3);
t43 = (-mrSges(4,1) * t65 + mrSges(4,2) * t62) * t58;
t81 = t61 * t17 + t84 * t18 + t39 * t95;
t82 = m(6) * (-pkin(4) * t68 + qJDD(3) * qJ(5) + 0.2e1 * qJD(5) * qJD(3) - t24 * t39 + t81) + qJD(3) * t36 + qJDD(3) * mrSges(6,3);
t9 = m(5) * t81 - qJDD(3) * mrSges(5,2) - qJD(3) * t35 + t31 * t87 + t39 * t86 + t82;
t6 = m(4) * (-t65 * g(3) - t88) - t44 * mrSges(4,3) + qJDD(3) * mrSges(4,1) - t43 * t90 + qJD(3) * t53 + t61 * t9 + t84 * t10;
t7 = m(4) * t79 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t45 - qJD(3) * t52 - t10 * t61 + t43 * t89 + t84 * t9;
t93 = t6 * t65 + t62 * t7;
t8 = m(3) * t78 + t57 * mrSges(3,1) - t56 * mrSges(3,2) - t96;
t3 = m(3) * t85 - mrSges(3,1) * t56 - mrSges(3,2) * t57 - t6 * t62 + t65 * t7;
t2 = m(2) * t77 - mrSges(2,1) * t69 - qJDD(1) * mrSges(2,2) + t3 * t66 - t63 * t8;
t1 = m(2) * t80 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t69 + t3 * t63 + t66 * t8;
t4 = [-m(1) * g(1) - t1 * t64 + t2 * t67, t2, t3, t7, t9, -t31 * mrSges(6,2) - t39 * t25 + t82; -m(1) * g(2) + t1 * t67 + t2 * t64, t1, t8, t6, t10, -t73; (-m(1) + t94) * g(3) + t93, g(3) * t94 + t93, -m(3) * g(3) + t93, t96, t72, -qJDD(3) * mrSges(6,1) + t32 * mrSges(6,2) - qJD(3) * t37 + t40 * t25 + t92;];
f_new = t4;
