% Calculate vector of cutting forces with Newton-Euler
% S5PRRPR5
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
% Datum: 2019-12-05 16:28
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRRPR5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR5_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR5_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR5_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR5_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR5_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:25:33
% EndTime: 2019-12-05 16:25:37
% DurationCPUTime: 1.44s
% Computational Cost: add. (15871->137), mult. (33578->193), div. (0->0), fcn. (23169->12), ass. (0->79)
t68 = sin(pkin(10));
t71 = cos(pkin(10));
t75 = sin(qJ(3));
t78 = cos(qJ(3));
t50 = (t68 * t75 - t71 * t78) * qJD(2);
t69 = sin(pkin(9));
t72 = cos(pkin(9));
t59 = t69 * g(1) - t72 * g(2);
t73 = cos(pkin(5));
t101 = t59 * t73;
t60 = -t72 * g(1) - t69 * g(2);
t67 = -g(3) + qJDD(1);
t70 = sin(pkin(5));
t76 = sin(qJ(2));
t79 = cos(qJ(2));
t105 = -t76 * t60 + (t67 * t70 + t101) * t79;
t96 = qJD(2) * qJD(3);
t92 = t78 * t96;
t57 = t75 * qJDD(2) + t92;
t58 = t78 * qJDD(2) - t75 * t96;
t98 = qJD(2) * t75;
t62 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t98;
t97 = qJD(2) * t78;
t63 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t97;
t81 = qJD(2) ^ 2;
t85 = -qJDD(2) * pkin(2) - t105;
t51 = (t68 * t78 + t71 * t75) * qJD(2);
t36 = t50 * pkin(4) - t51 * pkin(8);
t80 = qJD(3) ^ 2;
t103 = 2 * qJD(4);
t100 = t70 * t76;
t94 = t67 * t100 + t76 * t101 + t79 * t60;
t33 = -t81 * pkin(2) + qJDD(2) * pkin(7) + t94;
t47 = -t70 * t59 + t73 * t67;
t91 = -t75 * t33 + t78 * t47;
t22 = (-t57 + t92) * qJ(4) + (t75 * t78 * t81 + qJDD(3)) * pkin(3) + t91;
t61 = qJD(3) * pkin(3) - qJ(4) * t98;
t66 = t78 ^ 2;
t99 = t78 * t33 + t75 * t47;
t23 = -t66 * t81 * pkin(3) + t58 * qJ(4) - qJD(3) * t61 + t99;
t95 = -t50 * t103 + t68 * t22 + t71 * t23;
t18 = -t80 * pkin(4) + qJDD(3) * pkin(8) - t50 * t36 + t95;
t39 = -t68 * t57 + t71 * t58;
t40 = t71 * t57 + t68 * t58;
t82 = -t58 * pkin(3) + qJDD(4) + t61 * t98 + (-qJ(4) * t66 - pkin(7)) * t81 + t85;
t19 = (qJD(3) * t50 - t40) * pkin(8) + (qJD(3) * t51 - t39) * pkin(4) + t82;
t74 = sin(qJ(5));
t77 = cos(qJ(5));
t41 = t77 * qJD(3) - t74 * t51;
t27 = t41 * qJD(5) + t74 * qJDD(3) + t77 * t40;
t42 = t74 * qJD(3) + t77 * t51;
t28 = -t41 * mrSges(6,1) + t42 * mrSges(6,2);
t49 = qJD(5) + t50;
t30 = -t49 * mrSges(6,2) + t41 * mrSges(6,3);
t38 = qJDD(5) - t39;
t15 = m(6) * (-t74 * t18 + t77 * t19) - t27 * mrSges(6,3) + t38 * mrSges(6,1) - t42 * t28 + t49 * t30;
t26 = -t42 * qJD(5) + t77 * qJDD(3) - t74 * t40;
t31 = t49 * mrSges(6,1) - t42 * mrSges(6,3);
t16 = m(6) * (t77 * t18 + t74 * t19) + t26 * mrSges(6,3) - t38 * mrSges(6,2) + t41 * t28 - t49 * t31;
t45 = -qJD(3) * mrSges(5,2) - t50 * mrSges(5,3);
t46 = qJD(3) * mrSges(5,1) - t51 * mrSges(5,3);
t86 = -m(5) * t82 + t39 * mrSges(5,1) - t40 * mrSges(5,2) - t77 * t15 - t74 * t16 - t50 * t45 - t51 * t46;
t104 = (t75 * t62 - t78 * t63) * qJD(2) + m(4) * (-t81 * pkin(7) + t85) - t58 * mrSges(4,1) + t57 * mrSges(4,2) - t86;
t10 = m(3) * t105 + qJDD(2) * mrSges(3,1) - t81 * mrSges(3,2) - t104;
t102 = t10 * t79;
t35 = t50 * mrSges(5,1) + t51 * mrSges(5,2);
t11 = m(5) * t95 - qJDD(3) * mrSges(5,2) + t39 * mrSges(5,3) - qJD(3) * t46 - t74 * t15 + t77 * t16 - t50 * t35;
t90 = -t71 * t22 + t68 * t23;
t84 = m(6) * (-qJDD(3) * pkin(4) - t80 * pkin(8) + (t103 + t36) * t51 + t90) - t26 * mrSges(6,1) + t27 * mrSges(6,2) - t41 * t30 + t42 * t31;
t12 = m(5) * (-0.2e1 * qJD(4) * t51 - t90) - t40 * mrSges(5,3) + qJDD(3) * mrSges(5,1) - t51 * t35 + qJD(3) * t45 - t84;
t56 = (-mrSges(4,1) * t78 + mrSges(4,2) * t75) * qJD(2);
t7 = m(4) * t91 + qJDD(3) * mrSges(4,1) - t57 * mrSges(4,3) + qJD(3) * t63 + t68 * t11 + t71 * t12 - t56 * t98;
t8 = m(4) * t99 - qJDD(3) * mrSges(4,2) + t58 * mrSges(4,3) - qJD(3) * t62 + t71 * t11 - t68 * t12 + t56 * t97;
t4 = m(3) * t94 - t81 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t75 * t7 + t78 * t8;
t6 = m(3) * t47 + t78 * t7 + t75 * t8;
t93 = m(2) * t67 + t4 * t100 + t70 * t102 + t73 * t6;
t2 = m(2) * t60 - t76 * t10 + t79 * t4;
t1 = m(2) * t59 - t70 * t6 + (t4 * t76 + t102) * t73;
t3 = [-m(1) * g(1) - t69 * t1 + t72 * t2, t2, t4, t8, t11, t16; -m(1) * g(2) + t72 * t1 + t69 * t2, t1, t10, t7, t12, t15; -m(1) * g(3) + t93, t93, t6, t104, -t86, t84;];
f_new = t3;
