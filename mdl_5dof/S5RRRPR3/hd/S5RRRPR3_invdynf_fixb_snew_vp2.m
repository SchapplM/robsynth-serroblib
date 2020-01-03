% Calculate vector of cutting forces with Newton-Euler
% S5RRRPR3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2020-01-03 12:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRPR3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR3_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR3_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR3_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR3_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:08:56
% EndTime: 2020-01-03 12:08:58
% DurationCPUTime: 1.38s
% Computational Cost: add. (23161->141), mult. (31070->190), div. (0->0), fcn. (19687->10), ass. (0->76)
t64 = qJDD(1) + qJDD(2);
t71 = sin(qJ(3));
t75 = cos(qJ(3));
t66 = qJD(1) + qJD(2);
t91 = qJD(3) * t66;
t51 = t71 * t64 + t75 * t91;
t52 = t75 * t64 - t71 * t91;
t95 = t66 * t71;
t59 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t95;
t94 = t66 * t75;
t60 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t94;
t62 = t66 ^ 2;
t68 = sin(pkin(9));
t69 = cos(pkin(9));
t39 = -t68 * t51 + t69 * t52;
t40 = t69 * t51 + t68 * t52;
t47 = (-t68 * t71 + t69 * t75) * t66;
t41 = -qJD(3) * mrSges(5,2) + t47 * mrSges(5,3);
t48 = (t68 * t75 + t69 * t71) * t66;
t42 = qJD(3) * mrSges(5,1) - t48 * mrSges(5,3);
t58 = qJD(3) * pkin(3) - qJ(4) * t95;
t67 = t75 ^ 2;
t73 = sin(qJ(1));
t77 = cos(qJ(1));
t85 = -t77 * g(2) - t73 * g(3);
t56 = qJDD(1) * pkin(1) + t85;
t78 = qJD(1) ^ 2;
t88 = -t73 * g(2) + t77 * g(3);
t57 = -t78 * pkin(1) + t88;
t72 = sin(qJ(2));
t76 = cos(qJ(2));
t87 = t76 * t56 - t72 * t57;
t83 = -t64 * pkin(2) - t87;
t80 = -t52 * pkin(3) + qJDD(4) + t58 * t95 + (-qJ(4) * t67 - pkin(7)) * t62 + t83;
t70 = sin(qJ(5));
t74 = cos(qJ(5));
t32 = t70 * t47 + t74 * t48;
t18 = -t32 * qJD(5) + t74 * t39 - t70 * t40;
t31 = t74 * t47 - t70 * t48;
t19 = t31 * qJD(5) + t70 * t39 + t74 * t40;
t65 = qJD(3) + qJD(5);
t29 = -t65 * mrSges(6,2) + t31 * mrSges(6,3);
t30 = t65 * mrSges(6,1) - t32 * mrSges(6,3);
t43 = qJD(3) * pkin(4) - t48 * pkin(8);
t46 = t47 ^ 2;
t82 = t18 * mrSges(6,1) + t31 * t29 - m(6) * (-t39 * pkin(4) - t46 * pkin(8) + t48 * t43 + t80) - t19 * mrSges(6,2) - t32 * t30;
t81 = -m(5) * t80 + t39 * mrSges(5,1) - t40 * mrSges(5,2) + t47 * t41 - t48 * t42 + t82;
t99 = (t71 * t59 - t75 * t60) * t66 + m(4) * (-t62 * pkin(7) + t83) - t52 * mrSges(4,1) + t51 * mrSges(4,2) - t81;
t98 = -m(2) - m(3);
t50 = (-mrSges(4,1) * t75 + mrSges(4,2) * t71) * t66;
t92 = t72 * t56 + t76 * t57;
t37 = -t62 * pkin(2) + t64 * pkin(7) + t92;
t93 = t71 * t37;
t96 = pkin(3) * t62;
t25 = qJDD(3) * pkin(3) - t51 * qJ(4) - t93 + (qJ(4) * t91 + t71 * t96 - g(1)) * t75;
t89 = -t71 * g(1) + t75 * t37;
t26 = t52 * qJ(4) - qJD(3) * t58 - t67 * t96 + t89;
t86 = -0.2e1 * qJD(4) * t48 + t69 * t25 - t68 * t26;
t13 = (qJD(3) * t47 - t40) * pkin(8) + (t47 * t48 + qJDD(3)) * pkin(4) + t86;
t90 = 0.2e1 * qJD(4) * t47 + t68 * t25 + t69 * t26;
t14 = -t46 * pkin(4) + t39 * pkin(8) - qJD(3) * t43 + t90;
t21 = -t31 * mrSges(6,1) + t32 * mrSges(6,2);
t63 = qJDD(3) + qJDD(5);
t11 = m(6) * (t74 * t13 - t70 * t14) - t19 * mrSges(6,3) + t63 * mrSges(6,1) - t32 * t21 + t65 * t29;
t12 = m(6) * (t70 * t13 + t74 * t14) + t18 * mrSges(6,3) - t63 * mrSges(6,2) + t31 * t21 - t65 * t30;
t35 = -t47 * mrSges(5,1) + t48 * mrSges(5,2);
t8 = m(5) * t86 + qJDD(3) * mrSges(5,1) - t40 * mrSges(5,3) + qJD(3) * t41 + t74 * t11 + t70 * t12 - t48 * t35;
t9 = m(5) * t90 - qJDD(3) * mrSges(5,2) + t39 * mrSges(5,3) - qJD(3) * t42 - t70 * t11 + t74 * t12 + t47 * t35;
t6 = m(4) * (-t75 * g(1) - t93) - t51 * mrSges(4,3) + qJDD(3) * mrSges(4,1) - t50 * t95 + qJD(3) * t60 + t68 * t9 + t69 * t8;
t7 = m(4) * t89 - qJDD(3) * mrSges(4,2) + t52 * mrSges(4,3) - qJD(3) * t59 + t50 * t94 - t68 * t8 + t69 * t9;
t97 = t75 * t6 + t71 * t7;
t10 = m(3) * t87 + t64 * mrSges(3,1) - t62 * mrSges(3,2) - t99;
t3 = m(3) * t92 - t62 * mrSges(3,1) - t64 * mrSges(3,2) - t71 * t6 + t75 * t7;
t2 = m(2) * t85 + qJDD(1) * mrSges(2,1) - t78 * mrSges(2,2) + t76 * t10 + t72 * t3;
t1 = m(2) * t88 - t78 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t72 * t10 + t76 * t3;
t4 = [(-m(1) + t98) * g(1) + t97, t1, t3, t7, t9, t12; -m(1) * g(2) + t73 * t1 + t77 * t2, t2, t10, t6, t8, t11; -m(1) * g(3) - t77 * t1 + t73 * t2, t98 * g(1) + t97, -m(3) * g(1) + t97, t99, -t81, -t82;];
f_new = t4;
