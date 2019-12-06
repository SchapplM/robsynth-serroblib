% Calculate vector of cutting forces with Newton-Euler
% S5PRRPR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRRPR7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR7_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR7_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR7_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR7_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR7_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:35:14
% EndTime: 2019-12-05 16:35:18
% DurationCPUTime: 1.36s
% Computational Cost: add. (14460->137), mult. (29568->190), div. (0->0), fcn. (20141->12), ass. (0->81)
t65 = sin(pkin(9));
t68 = cos(pkin(9));
t56 = -t68 * g(1) - t65 * g(2);
t63 = -g(3) + qJDD(1);
t66 = sin(pkin(5));
t72 = sin(qJ(2));
t75 = cos(qJ(2));
t55 = t65 * g(1) - t68 * g(2);
t69 = cos(pkin(5));
t97 = t55 * t69;
t101 = -t72 * t56 + (t63 * t66 + t97) * t75;
t64 = sin(pkin(10));
t67 = cos(pkin(10));
t71 = sin(qJ(3));
t92 = qJD(2) * t71;
t47 = t67 * qJD(3) - t64 * t92;
t48 = t64 * qJD(3) + t67 * t92;
t33 = -t47 * pkin(4) - t48 * pkin(8);
t74 = cos(qJ(3));
t90 = qJD(2) * qJD(3);
t86 = t71 * t90;
t54 = t74 * qJDD(2) - t86;
t51 = (-pkin(3) * t74 - qJ(4) * t71) * qJD(2);
t76 = qJD(3) ^ 2;
t91 = qJD(2) * t74;
t77 = qJD(2) ^ 2;
t95 = t66 * t72;
t88 = t75 * t56 + t63 * t95 + t72 * t97;
t29 = -t77 * pkin(2) + qJDD(2) * pkin(7) + t88;
t42 = -t66 * t55 + t69 * t63;
t93 = t74 * t29 + t71 * t42;
t20 = -t76 * pkin(3) + qJDD(3) * qJ(4) + t51 * t91 + t93;
t28 = -qJDD(2) * pkin(2) - t77 * pkin(7) - t101;
t85 = t74 * t90;
t53 = t71 * qJDD(2) + t85;
t22 = (-t53 - t85) * qJ(4) + (-t54 + t86) * pkin(3) + t28;
t99 = 2 * qJD(4);
t89 = t67 * t20 + t64 * t22 + t47 * t99;
t96 = t74 ^ 2 * t77;
t16 = -pkin(4) * t96 - t54 * pkin(8) + t47 * t33 + t89;
t40 = t67 * qJDD(3) - t64 * t53;
t41 = t64 * qJDD(3) + t67 * t53;
t26 = t71 * t29;
t81 = -qJDD(3) * pkin(3) - t76 * qJ(4) + t51 * t92 + qJDD(4) + t26;
t17 = -t40 * pkin(4) - t41 * pkin(8) + (-t42 + (-pkin(4) * t48 + pkin(8) * t47) * qJD(2)) * t74 + t81;
t70 = sin(qJ(5));
t73 = cos(qJ(5));
t34 = -t70 * t48 - t73 * t91;
t24 = t34 * qJD(5) + t73 * t41 - t70 * t54;
t35 = t73 * t48 - t70 * t91;
t25 = -t34 * mrSges(6,1) + t35 * mrSges(6,2);
t46 = qJD(5) - t47;
t30 = -t46 * mrSges(6,2) + t34 * mrSges(6,3);
t37 = qJDD(5) - t40;
t13 = m(6) * (-t70 * t16 + t73 * t17) - t24 * mrSges(6,3) + t37 * mrSges(6,1) - t35 * t25 + t46 * t30;
t23 = -t35 * qJD(5) - t70 * t41 - t73 * t54;
t31 = t46 * mrSges(6,1) - t35 * mrSges(6,3);
t14 = m(6) * (t73 * t16 + t70 * t17) + t23 * mrSges(6,3) - t37 * mrSges(6,2) + t34 * t25 - t46 * t31;
t32 = -t47 * mrSges(5,1) + t48 * mrSges(5,2);
t39 = -mrSges(5,1) * t91 - t48 * mrSges(5,3);
t11 = m(5) * t89 + t54 * mrSges(5,2) + t40 * mrSges(5,3) - t70 * t13 + t73 * t14 + t47 * t32 + t39 * t91;
t38 = mrSges(5,2) * t91 + t47 * mrSges(5,3);
t84 = t64 * t20 - t67 * t22;
t79 = m(6) * (-pkin(8) * t96 + t54 * pkin(4) + (t99 + t33) * t48 + t84) - t23 * mrSges(6,1) + t24 * mrSges(6,2) - t34 * t30 + t35 * t31;
t12 = m(5) * (-0.2e1 * qJD(4) * t48 - t84) - t41 * mrSges(5,3) - t54 * mrSges(5,1) - t48 * t32 - t38 * t91 - t79;
t57 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t92;
t58 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t91;
t100 = m(4) * t28 - t54 * mrSges(4,1) + t53 * mrSges(4,2) + t64 * t11 + t67 * t12 + (t71 * t57 - t74 * t58) * qJD(2);
t8 = m(3) * t101 + qJDD(2) * mrSges(3,1) - t77 * mrSges(3,2) - t100;
t98 = t75 * t8;
t94 = t74 * t42;
t52 = (-mrSges(4,1) * t74 + mrSges(4,2) * t71) * qJD(2);
t78 = m(5) * (t81 - t94) - t40 * mrSges(5,1) + t41 * mrSges(5,2) + t73 * t13 + t70 * t14 - t47 * t38 + t48 * t39;
t10 = m(4) * (-t26 + t94) - t53 * mrSges(4,3) + qJDD(3) * mrSges(4,1) - t52 * t92 + qJD(3) * t58 - t78;
t9 = m(4) * t93 - qJDD(3) * mrSges(4,2) + t54 * mrSges(4,3) - qJD(3) * t57 + t67 * t11 - t64 * t12 + t52 * t91;
t4 = m(3) * t88 - t77 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t71 * t10 + t74 * t9;
t6 = m(3) * t42 + t74 * t10 + t71 * t9;
t87 = m(2) * t63 + t4 * t95 + t69 * t6 + t66 * t98;
t2 = m(2) * t56 + t75 * t4 - t72 * t8;
t1 = m(2) * t55 - t66 * t6 + (t4 * t72 + t98) * t69;
t3 = [-m(1) * g(1) - t65 * t1 + t68 * t2, t2, t4, t9, t11, t14; -m(1) * g(2) + t68 * t1 + t65 * t2, t1, t8, t10, t12, t13; -m(1) * g(3) + t87, t87, t6, t100, t78, t79;];
f_new = t3;
