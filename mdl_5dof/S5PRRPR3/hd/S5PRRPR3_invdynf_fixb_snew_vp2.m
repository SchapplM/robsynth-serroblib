% Calculate vector of cutting forces with Newton-Euler
% S5PRRPR3
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRRPR3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR3_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR3_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR3_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR3_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR3_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:19:09
% EndTime: 2019-12-05 16:19:12
% DurationCPUTime: 1.25s
% Computational Cost: add. (13369->131), mult. (29459->182), div. (0->0), fcn. (19687->10), ass. (0->72)
t74 = sin(qJ(3));
t77 = cos(qJ(3));
t93 = qJD(2) * qJD(3);
t90 = t77 * t93;
t54 = t74 * qJDD(2) + t90;
t55 = t77 * qJDD(2) - t74 * t93;
t95 = qJD(2) * t74;
t59 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t95;
t94 = qJD(2) * t77;
t60 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t94;
t79 = qJD(2) ^ 2;
t69 = sin(pkin(9));
t71 = cos(pkin(9));
t39 = -t69 * t54 + t71 * t55;
t40 = t71 * t54 + t69 * t55;
t47 = (-t69 * t74 + t71 * t77) * qJD(2);
t41 = -qJD(3) * mrSges(5,2) + t47 * mrSges(5,3);
t48 = (t69 * t77 + t71 * t74) * qJD(2);
t42 = qJD(3) * mrSges(5,1) - t48 * mrSges(5,3);
t58 = qJD(3) * pkin(3) - qJ(4) * t95;
t67 = t77 ^ 2;
t70 = sin(pkin(8));
t72 = cos(pkin(8));
t56 = t70 * g(1) - t72 * g(2);
t57 = -t72 * g(1) - t70 * g(2);
t75 = sin(qJ(2));
t78 = cos(qJ(2));
t87 = t78 * t56 - t75 * t57;
t84 = -qJDD(2) * pkin(2) - t87;
t81 = -t55 * pkin(3) + qJDD(4) + t58 * t95 + (-qJ(4) * t67 - pkin(6)) * t79 + t84;
t73 = sin(qJ(5));
t76 = cos(qJ(5));
t32 = t73 * t47 + t76 * t48;
t18 = -t32 * qJD(5) + t76 * t39 - t73 * t40;
t31 = t76 * t47 - t73 * t48;
t19 = t31 * qJD(5) + t73 * t39 + t76 * t40;
t66 = qJD(3) + qJD(5);
t29 = -t66 * mrSges(6,2) + t31 * mrSges(6,3);
t30 = t66 * mrSges(6,1) - t32 * mrSges(6,3);
t43 = qJD(3) * pkin(4) - t48 * pkin(7);
t46 = t47 ^ 2;
t83 = t18 * mrSges(6,1) + t31 * t29 - m(6) * (-t39 * pkin(4) - t46 * pkin(7) + t48 * t43 + t81) - t19 * mrSges(6,2) - t32 * t30;
t82 = -m(5) * t81 + t39 * mrSges(5,1) - t40 * mrSges(5,2) + t47 * t41 - t48 * t42 + t83;
t98 = (t74 * t59 - t77 * t60) * qJD(2) + m(4) * (-t79 * pkin(6) + t84) - t55 * mrSges(4,1) + t54 * mrSges(4,2) - t82;
t96 = t75 * t56 + t78 * t57;
t38 = -t79 * pkin(2) + qJDD(2) * pkin(6) + t96;
t68 = -g(3) + qJDD(1);
t97 = t77 * t38 + t74 * t68;
t53 = (-mrSges(4,1) * t77 + mrSges(4,2) * t74) * qJD(2);
t88 = -t74 * t38 + t77 * t68;
t25 = (-t54 + t90) * qJ(4) + (t74 * t77 * t79 + qJDD(3)) * pkin(3) + t88;
t26 = -t67 * t79 * pkin(3) + t55 * qJ(4) - qJD(3) * t58 + t97;
t86 = -0.2e1 * qJD(4) * t48 + t71 * t25 - t69 * t26;
t13 = (qJD(3) * t47 - t40) * pkin(7) + (t47 * t48 + qJDD(3)) * pkin(4) + t86;
t91 = 0.2e1 * qJD(4) * t47 + t69 * t25 + t71 * t26;
t14 = -t46 * pkin(4) + t39 * pkin(7) - qJD(3) * t43 + t91;
t21 = -t31 * mrSges(6,1) + t32 * mrSges(6,2);
t65 = qJDD(3) + qJDD(5);
t11 = m(6) * (t76 * t13 - t73 * t14) - t19 * mrSges(6,3) + t65 * mrSges(6,1) - t32 * t21 + t66 * t29;
t12 = m(6) * (t73 * t13 + t76 * t14) + t18 * mrSges(6,3) - t65 * mrSges(6,2) + t31 * t21 - t66 * t30;
t34 = -t47 * mrSges(5,1) + t48 * mrSges(5,2);
t8 = m(5) * t86 + qJDD(3) * mrSges(5,1) - t40 * mrSges(5,3) + qJD(3) * t41 + t76 * t11 + t73 * t12 - t48 * t34;
t9 = m(5) * t91 - qJDD(3) * mrSges(5,2) + t39 * mrSges(5,3) - qJD(3) * t42 - t73 * t11 + t76 * t12 + t47 * t34;
t6 = m(4) * t88 + qJDD(3) * mrSges(4,1) - t54 * mrSges(4,3) + qJD(3) * t60 - t53 * t95 + t69 * t9 + t71 * t8;
t7 = m(4) * t97 - qJDD(3) * mrSges(4,2) + t55 * mrSges(4,3) - qJD(3) * t59 + t53 * t94 - t69 * t8 + t71 * t9;
t92 = m(3) * t68 + t77 * t6 + t74 * t7;
t89 = m(2) * t68 + t92;
t10 = m(3) * t87 + qJDD(2) * mrSges(3,1) - t79 * mrSges(3,2) - t98;
t3 = m(3) * t96 - t79 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t74 * t6 + t77 * t7;
t2 = m(2) * t57 - t75 * t10 + t78 * t3;
t1 = m(2) * t56 + t78 * t10 + t75 * t3;
t4 = [-m(1) * g(1) - t70 * t1 + t72 * t2, t2, t3, t7, t9, t12; -m(1) * g(2) + t72 * t1 + t70 * t2, t1, t10, t6, t8, t11; -m(1) * g(3) + t89, t89, t92, t98, -t82, -t83;];
f_new = t4;
