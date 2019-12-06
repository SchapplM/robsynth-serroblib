% Calculate vector of cutting forces with Newton-Euler
% S5PRRPR6
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
% Datum: 2019-12-05 16:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRRPR6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR6_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR6_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR6_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR6_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR6_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:30:21
% EndTime: 2019-12-05 16:30:25
% DurationCPUTime: 1.52s
% Computational Cost: add. (17476->135), mult. (35758->188), div. (0->0), fcn. (24513->12), ass. (0->79)
t70 = sin(pkin(9));
t73 = cos(pkin(9));
t61 = -t73 * g(1) - t70 * g(2);
t68 = -g(3) + qJDD(1);
t71 = sin(pkin(5));
t77 = sin(qJ(2));
t80 = cos(qJ(2));
t60 = t70 * g(1) - t73 * g(2);
t74 = cos(pkin(5));
t99 = t60 * t74;
t102 = -t77 * t61 + (t68 * t71 + t99) * t80;
t76 = sin(qJ(3));
t79 = cos(qJ(3));
t94 = qJD(2) * qJD(3);
t90 = t79 * t94;
t58 = t76 * qJDD(2) + t90;
t69 = sin(pkin(10));
t72 = cos(pkin(10));
t43 = t69 * qJDD(3) + t72 * t58;
t96 = qJD(2) * t76;
t51 = t72 * qJD(3) - t69 * t96;
t52 = t69 * qJD(3) + t72 * t96;
t66 = t76 * t94;
t59 = t79 * qJDD(2) - t66;
t56 = (-pkin(3) * t79 - qJ(4) * t76) * qJD(2);
t81 = qJD(3) ^ 2;
t95 = t79 * qJD(2);
t82 = qJD(2) ^ 2;
t98 = t71 * t77;
t92 = t80 * t61 + t68 * t98 + t77 * t99;
t33 = -t82 * pkin(2) + qJDD(2) * pkin(7) + t92;
t45 = -t71 * t60 + t74 * t68;
t97 = t79 * t33 + t76 * t45;
t21 = -t81 * pkin(3) + qJDD(3) * qJ(4) + t56 * t95 + t97;
t32 = -qJDD(2) * pkin(2) - t82 * pkin(7) - t102;
t25 = (-t58 - t90) * qJ(4) + (-t59 + t66) * pkin(3) + t32;
t88 = -0.2e1 * qJD(4) * t52 - t69 * t21 + t72 * t25;
t15 = (-t51 * t95 - t43) * pkin(8) + (t51 * t52 - t59) * pkin(4) + t88;
t42 = t72 * qJDD(3) - t69 * t58;
t44 = -pkin(4) * t95 - t52 * pkin(8);
t50 = t51 ^ 2;
t93 = 0.2e1 * qJD(4) * t51 + t72 * t21 + t69 * t25;
t16 = -t50 * pkin(4) + t42 * pkin(8) + t44 * t95 + t93;
t75 = sin(qJ(5));
t78 = cos(qJ(5));
t36 = t78 * t51 - t75 * t52;
t27 = t36 * qJD(5) + t75 * t42 + t78 * t43;
t37 = t75 * t51 + t78 * t52;
t29 = -t36 * mrSges(6,1) + t37 * mrSges(6,2);
t65 = qJD(5) - t95;
t34 = -t65 * mrSges(6,2) + t36 * mrSges(6,3);
t53 = qJDD(5) - t59;
t13 = m(6) * (t78 * t15 - t75 * t16) - t27 * mrSges(6,3) + t53 * mrSges(6,1) - t37 * t29 + t65 * t34;
t26 = -t37 * qJD(5) + t78 * t42 - t75 * t43;
t35 = t65 * mrSges(6,1) - t37 * mrSges(6,3);
t14 = m(6) * (t75 * t15 + t78 * t16) + t26 * mrSges(6,3) - t53 * mrSges(6,2) + t36 * t29 - t65 * t35;
t38 = -t51 * mrSges(5,1) + t52 * mrSges(5,2);
t40 = mrSges(5,2) * t95 + t51 * mrSges(5,3);
t10 = m(5) * t88 - t59 * mrSges(5,1) - t43 * mrSges(5,3) + t78 * t13 + t75 * t14 - t52 * t38 - t40 * t95;
t41 = -mrSges(5,1) * t95 - t52 * mrSges(5,3);
t11 = m(5) * t93 + t59 * mrSges(5,2) + t42 * mrSges(5,3) - t75 * t13 + t78 * t14 + t51 * t38 + t41 * t95;
t62 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t96;
t63 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t95;
t101 = m(4) * t32 - t59 * mrSges(4,1) + t58 * mrSges(4,2) + t72 * t10 + t69 * t11 + (t76 * t62 - t79 * t63) * qJD(2);
t8 = m(3) * t102 + qJDD(2) * mrSges(3,1) - t82 * mrSges(3,2) - t101;
t100 = t8 * t80;
t57 = (-mrSges(4,1) * t79 + mrSges(4,2) * t76) * qJD(2);
t89 = -t76 * t33 + t79 * t45;
t20 = -qJDD(3) * pkin(3) - t81 * qJ(4) + t56 * t96 + qJDD(4) - t89;
t85 = t26 * mrSges(6,1) + t36 * t34 - m(6) * (-t42 * pkin(4) - t50 * pkin(8) + t52 * t44 + t20) - t27 * mrSges(6,2) - t37 * t35;
t83 = m(5) * t20 - t42 * mrSges(5,1) + t43 * mrSges(5,2) - t51 * t40 + t52 * t41 - t85;
t12 = m(4) * t89 + qJDD(3) * mrSges(4,1) - t58 * mrSges(4,3) + qJD(3) * t63 - t57 * t96 - t83;
t9 = m(4) * t97 - qJDD(3) * mrSges(4,2) + t59 * mrSges(4,3) - qJD(3) * t62 - t69 * t10 + t72 * t11 + t57 * t95;
t4 = m(3) * t92 - t82 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t76 * t12 + t79 * t9;
t6 = m(3) * t45 + t79 * t12 + t76 * t9;
t91 = m(2) * t68 + t71 * t100 + t4 * t98 + t74 * t6;
t2 = m(2) * t61 + t80 * t4 - t77 * t8;
t1 = m(2) * t60 - t71 * t6 + (t4 * t77 + t100) * t74;
t3 = [-m(1) * g(1) - t70 * t1 + t73 * t2, t2, t4, t9, t11, t14; -m(1) * g(2) + t73 * t1 + t70 * t2, t1, t8, t12, t10, t13; -m(1) * g(3) + t91, t91, t6, t101, t83, -t85;];
f_new = t3;
