% Calculate vector of cutting forces with Newton-Euler
% S5PRRRR9
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRRRR9_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR9_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR9_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR9_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR9_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR9_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR9_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR9_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR9_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:19:15
% EndTime: 2019-12-05 17:19:20
% DurationCPUTime: 1.57s
% Computational Cost: add. (19149->136), mult. (36860->186), div. (0->0), fcn. (25545->12), ass. (0->81)
t69 = sin(pkin(10));
t71 = cos(pkin(10));
t60 = -t71 * g(1) - t69 * g(2);
t68 = -g(3) + qJDD(1);
t70 = sin(pkin(5));
t76 = sin(qJ(2));
t80 = cos(qJ(2));
t59 = t69 * g(1) - t71 * g(2);
t72 = cos(pkin(5));
t99 = t59 * t72;
t102 = -t76 * t60 + (t68 * t70 + t99) * t80;
t74 = sin(qJ(4));
t78 = cos(qJ(4));
t75 = sin(qJ(3));
t95 = qJD(2) * t75;
t53 = t78 * qJD(3) - t74 * t95;
t79 = cos(qJ(3));
t93 = qJD(2) * qJD(3);
t90 = t79 * t93;
t57 = t75 * qJDD(2) + t90;
t37 = t53 * qJD(4) + t74 * qJDD(3) + t78 * t57;
t66 = t75 * t93;
t58 = t79 * qJDD(2) - t66;
t50 = qJDD(4) - t58;
t54 = t74 * qJD(3) + t78 * t95;
t94 = t79 * qJD(2);
t65 = qJD(4) - t94;
t56 = (-pkin(3) * t79 - pkin(8) * t75) * qJD(2);
t81 = qJD(3) ^ 2;
t82 = qJD(2) ^ 2;
t98 = t70 * t76;
t92 = t80 * t60 + t68 * t98 + t76 * t99;
t33 = -t82 * pkin(2) + qJDD(2) * pkin(7) + t92;
t45 = -t70 * t59 + t72 * t68;
t96 = t79 * t33 + t75 * t45;
t24 = -t81 * pkin(3) + qJDD(3) * pkin(8) + t56 * t94 + t96;
t32 = -qJDD(2) * pkin(2) - t82 * pkin(7) - t102;
t27 = (-t57 - t90) * pkin(8) + (-t58 + t66) * pkin(3) + t32;
t89 = -t74 * t24 + t78 * t27;
t15 = (t53 * t65 - t37) * pkin(9) + (t53 * t54 + t50) * pkin(4) + t89;
t36 = -t54 * qJD(4) + t78 * qJDD(3) - t74 * t57;
t44 = t65 * pkin(4) - t54 * pkin(9);
t49 = t53 ^ 2;
t97 = t78 * t24 + t74 * t27;
t16 = -t49 * pkin(4) + t36 * pkin(9) - t65 * t44 + t97;
t73 = sin(qJ(5));
t77 = cos(qJ(5));
t38 = t77 * t53 - t73 * t54;
t21 = t38 * qJD(5) + t73 * t36 + t77 * t37;
t39 = t73 * t53 + t77 * t54;
t29 = -t38 * mrSges(6,1) + t39 * mrSges(6,2);
t64 = qJD(5) + t65;
t34 = -t64 * mrSges(6,2) + t38 * mrSges(6,3);
t48 = qJDD(5) + t50;
t13 = m(6) * (t77 * t15 - t73 * t16) - t21 * mrSges(6,3) + t48 * mrSges(6,1) - t39 * t29 + t64 * t34;
t20 = -t39 * qJD(5) + t77 * t36 - t73 * t37;
t35 = t64 * mrSges(6,1) - t39 * mrSges(6,3);
t14 = m(6) * (t73 * t15 + t77 * t16) + t20 * mrSges(6,3) - t48 * mrSges(6,2) + t38 * t29 - t64 * t35;
t40 = -t53 * mrSges(5,1) + t54 * mrSges(5,2);
t42 = -t65 * mrSges(5,2) + t53 * mrSges(5,3);
t10 = m(5) * t89 + t50 * mrSges(5,1) - t37 * mrSges(5,3) + t77 * t13 + t73 * t14 - t54 * t40 + t65 * t42;
t43 = t65 * mrSges(5,1) - t54 * mrSges(5,3);
t11 = m(5) * t97 - t50 * mrSges(5,2) + t36 * mrSges(5,3) - t73 * t13 + t77 * t14 + t53 * t40 - t65 * t43;
t61 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t95;
t62 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t94;
t101 = m(4) * t32 - t58 * mrSges(4,1) + t57 * mrSges(4,2) + t78 * t10 + t74 * t11 + (t75 * t61 - t79 * t62) * qJD(2);
t8 = m(3) * t102 + qJDD(2) * mrSges(3,1) - t82 * mrSges(3,2) - t101;
t100 = t8 * t80;
t55 = (-mrSges(4,1) * t79 + mrSges(4,2) * t75) * qJD(2);
t88 = -t75 * t33 + t79 * t45;
t23 = -qJDD(3) * pkin(3) - t81 * pkin(8) + t56 * t95 - t88;
t85 = t20 * mrSges(6,1) + t38 * t34 - m(6) * (-t36 * pkin(4) - t49 * pkin(9) + t54 * t44 + t23) - t21 * mrSges(6,2) - t39 * t35;
t83 = m(5) * t23 - t36 * mrSges(5,1) + t37 * mrSges(5,2) - t53 * t42 + t54 * t43 - t85;
t12 = m(4) * t88 + qJDD(3) * mrSges(4,1) - t57 * mrSges(4,3) + qJD(3) * t62 - t55 * t95 - t83;
t9 = m(4) * t96 - qJDD(3) * mrSges(4,2) + t58 * mrSges(4,3) - qJD(3) * t61 - t74 * t10 + t78 * t11 + t55 * t94;
t4 = m(3) * t92 - t82 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t75 * t12 + t79 * t9;
t6 = m(3) * t45 + t79 * t12 + t75 * t9;
t91 = m(2) * t68 + t70 * t100 + t4 * t98 + t72 * t6;
t2 = m(2) * t60 + t80 * t4 - t76 * t8;
t1 = m(2) * t59 - t70 * t6 + (t4 * t76 + t100) * t72;
t3 = [-m(1) * g(1) - t69 * t1 + t71 * t2, t2, t4, t9, t11, t14; -m(1) * g(2) + t71 * t1 + t69 * t2, t1, t8, t12, t10, t13; -m(1) * g(3) + t91, t91, t6, t101, t83, -t85;];
f_new = t3;
