% Calculate vector of cutting forces with Newton-Euler
% S5RPRPR8
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-31 18:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRPR8_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR8_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR8_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR8_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR8_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR8_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR8_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR8_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR8_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:21:23
% EndTime: 2019-12-31 18:21:25
% DurationCPUTime: 1.12s
% Computational Cost: add. (12443->136), mult. (25485->183), div. (0->0), fcn. (15701->10), ass. (0->74)
t77 = qJD(1) ^ 2;
t72 = sin(qJ(1));
t75 = cos(qJ(1));
t87 = t72 * g(1) - t75 * g(2);
t52 = qJDD(1) * pkin(1) + t87;
t82 = -t75 * g(1) - t72 * g(2);
t55 = -t77 * pkin(1) + t82;
t67 = sin(pkin(8));
t69 = cos(pkin(8));
t84 = t69 * t52 - t67 * t55;
t32 = -qJDD(1) * pkin(2) - t77 * pkin(6) - t84;
t71 = sin(qJ(3));
t74 = cos(qJ(3));
t90 = qJD(1) * qJD(3);
t86 = t74 * t90;
t56 = t71 * qJDD(1) + t86;
t62 = t71 * t90;
t57 = t74 * qJDD(1) - t62;
t92 = qJD(1) * t71;
t58 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t92;
t91 = t74 * qJD(1);
t59 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t91;
t66 = sin(pkin(9));
t68 = cos(pkin(9));
t40 = t66 * qJDD(3) + t68 * t56;
t49 = t68 * qJD(3) - t66 * t92;
t50 = t66 * qJD(3) + t68 * t92;
t23 = (-t56 - t86) * qJ(4) + (-t57 + t62) * pkin(3) + t32;
t53 = (-pkin(3) * t74 - qJ(4) * t71) * qJD(1);
t76 = qJD(3) ^ 2;
t93 = t67 * t52 + t69 * t55;
t33 = -t77 * pkin(2) + qJDD(1) * pkin(6) + t93;
t65 = -g(3) + qJDD(2);
t94 = t74 * t33 + t71 * t65;
t27 = -t76 * pkin(3) + qJDD(3) * qJ(4) + t53 * t91 + t94;
t83 = -0.2e1 * qJD(4) * t50 + t68 * t23 - t66 * t27;
t13 = (-t49 * t91 - t40) * pkin(7) + (t49 * t50 - t57) * pkin(4) + t83;
t39 = t68 * qJDD(3) - t66 * t56;
t41 = -pkin(4) * t91 - t50 * pkin(7);
t48 = t49 ^ 2;
t88 = 0.2e1 * qJD(4) * t49 + t66 * t23 + t68 * t27;
t14 = -t48 * pkin(4) + t39 * pkin(7) + t41 * t91 + t88;
t70 = sin(qJ(5));
t73 = cos(qJ(5));
t34 = t73 * t49 - t70 * t50;
t19 = t34 * qJD(5) + t70 * t39 + t73 * t40;
t35 = t70 * t49 + t73 * t50;
t26 = -t34 * mrSges(6,1) + t35 * mrSges(6,2);
t60 = qJD(5) - t91;
t28 = -t60 * mrSges(6,2) + t34 * mrSges(6,3);
t51 = qJDD(5) - t57;
t11 = m(6) * (t73 * t13 - t70 * t14) - t19 * mrSges(6,3) + t51 * mrSges(6,1) - t35 * t26 + t60 * t28;
t18 = -t35 * qJD(5) + t73 * t39 - t70 * t40;
t29 = t60 * mrSges(6,1) - t35 * mrSges(6,3);
t12 = m(6) * (t70 * t13 + t73 * t14) + t18 * mrSges(6,3) - t51 * mrSges(6,2) + t34 * t26 - t60 * t29;
t36 = -t49 * mrSges(5,1) + t50 * mrSges(5,2);
t37 = mrSges(5,2) * t91 + t49 * mrSges(5,3);
t7 = m(5) * t83 - t57 * mrSges(5,1) - t40 * mrSges(5,3) + t73 * t11 + t70 * t12 - t50 * t36 - t37 * t91;
t38 = -mrSges(5,1) * t91 - t50 * mrSges(5,3);
t8 = m(5) * t88 + t57 * mrSges(5,2) + t39 * mrSges(5,3) - t70 * t11 + t73 * t12 + t49 * t36 + t38 * t91;
t95 = m(4) * t32 - t57 * mrSges(4,1) + t56 * mrSges(4,2) + t66 * t8 + t68 * t7 + (t71 * t58 - t74 * t59) * qJD(1);
t54 = (-mrSges(4,1) * t74 + mrSges(4,2) * t71) * qJD(1);
t85 = -t71 * t33 + t74 * t65;
t25 = -qJDD(3) * pkin(3) - t76 * qJ(4) + t53 * t92 + qJDD(4) - t85;
t80 = t18 * mrSges(6,1) + t34 * t28 - m(6) * (-t39 * pkin(4) - t48 * pkin(7) + t50 * t41 + t25) - t19 * mrSges(6,2) - t35 * t29;
t78 = m(5) * t25 - t39 * mrSges(5,1) + t40 * mrSges(5,2) - t49 * t37 + t50 * t38 - t80;
t10 = m(4) * t85 + qJDD(3) * mrSges(4,1) - t56 * mrSges(4,3) + qJD(3) * t59 - t54 * t92 - t78;
t6 = m(4) * t94 - qJDD(3) * mrSges(4,2) + t57 * mrSges(4,3) - qJD(3) * t58 + t54 * t91 - t66 * t7 + t68 * t8;
t89 = m(3) * t65 + t74 * t10 + t71 * t6;
t4 = m(3) * t84 + qJDD(1) * mrSges(3,1) - t77 * mrSges(3,2) - t95;
t3 = m(3) * t93 - t77 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t71 * t10 + t74 * t6;
t2 = m(2) * t82 - t77 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t69 * t3 - t67 * t4;
t1 = m(2) * t87 + qJDD(1) * mrSges(2,1) - t77 * mrSges(2,2) + t67 * t3 + t69 * t4;
t5 = [-m(1) * g(1) - t72 * t1 + t75 * t2, t2, t3, t6, t8, t12; -m(1) * g(2) + t75 * t1 + t72 * t2, t1, t4, t10, t7, t11; (-m(1) - m(2)) * g(3) + t89, -m(2) * g(3) + t89, t89, t95, t78, -t80;];
f_new = t5;
