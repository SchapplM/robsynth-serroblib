% Calculate vector of cutting forces with Newton-Euler
% S5RPRPR12
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
% Datum: 2019-12-31 18:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRPR12_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR12_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR12_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR12_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR12_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR12_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR12_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR12_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR12_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:29:30
% EndTime: 2019-12-31 18:29:34
% DurationCPUTime: 1.88s
% Computational Cost: add. (19302->152), mult. (47477->202), div. (0->0), fcn. (34378->10), ass. (0->84)
t111 = cos(qJ(3));
t77 = sin(pkin(8));
t79 = cos(pkin(8));
t81 = sin(qJ(3));
t115 = -t79 * t111 + t77 * t81;
t86 = qJD(1) ^ 2;
t75 = t79 ^ 2;
t108 = t77 ^ 2 + t75;
t114 = t108 * mrSges(3,3);
t104 = qJD(1) * qJD(2);
t100 = -t79 * g(3) - 0.2e1 * t77 * t104;
t65 = t115 * qJD(1);
t92 = t111 * t77 + t79 * t81;
t66 = t92 * qJD(1);
t49 = t65 * mrSges(4,1) + t66 * mrSges(4,2);
t106 = t65 * qJD(3);
t54 = t92 * qJDD(1) - t106;
t60 = -qJD(3) * mrSges(4,2) - t65 * mrSges(4,3);
t48 = t65 * pkin(3) - t66 * qJ(4);
t85 = qJD(3) ^ 2;
t107 = pkin(6) * qJDD(1);
t112 = pkin(2) * t86;
t82 = sin(qJ(1));
t84 = cos(qJ(1));
t97 = -t84 * g(1) - t82 * g(2);
t67 = -t86 * pkin(1) + qJDD(1) * qJ(2) + t97;
t39 = (t79 * t112 - t107 - t67) * t77 + t100;
t99 = -t77 * g(3) + (0.2e1 * t104 + t67) * t79;
t45 = t79 * t107 - t75 * t112 + t99;
t95 = t111 * t39 - t81 * t45;
t21 = -qJDD(3) * pkin(3) - t85 * qJ(4) + t66 * t48 + qJDD(4) - t95;
t76 = sin(pkin(9));
t78 = cos(pkin(9));
t58 = t78 * qJD(3) - t76 * t66;
t40 = -t65 * mrSges(5,2) + t58 * mrSges(5,3);
t59 = t76 * qJD(3) + t78 * t66;
t41 = t65 * mrSges(5,1) - t59 * mrSges(5,3);
t43 = t78 * qJDD(3) - t76 * t54;
t44 = t76 * qJDD(3) + t78 * t54;
t80 = sin(qJ(5));
t83 = cos(qJ(5));
t33 = t80 * t58 + t83 * t59;
t22 = -t33 * qJD(5) + t83 * t43 - t80 * t44;
t32 = t83 * t58 - t80 * t59;
t23 = t32 * qJD(5) + t80 * t43 + t83 * t44;
t63 = qJD(5) + t65;
t30 = -t63 * mrSges(6,2) + t32 * mrSges(6,3);
t31 = t63 * mrSges(6,1) - t33 * mrSges(6,3);
t42 = t65 * pkin(4) - t59 * pkin(7);
t57 = t58 ^ 2;
t91 = t22 * mrSges(6,1) + t32 * t30 - m(6) * (-t43 * pkin(4) - t57 * pkin(7) + t59 * t42 + t21) - t23 * mrSges(6,2) - t33 * t31;
t87 = m(5) * t21 - t43 * mrSges(5,1) + t44 * mrSges(5,2) - t58 * t40 + t59 * t41 - t91;
t12 = m(4) * t95 + qJDD(3) * mrSges(4,1) - t54 * mrSges(4,3) + qJD(3) * t60 - t66 * t49 - t87;
t105 = t66 * qJD(3);
t53 = t115 * qJDD(1) + t105;
t109 = t111 * t45 + t81 * t39;
t24 = -t85 * pkin(3) + qJDD(3) * qJ(4) - t65 * t48 + t109;
t101 = t82 * g(1) - t84 * g(2);
t96 = qJDD(2) - t101;
t88 = (-pkin(2) * t79 - pkin(1)) * qJDD(1) + (-t108 * pkin(6) - qJ(2)) * t86 + t96;
t27 = (-t54 + t106) * qJ(4) + (t53 + t105) * pkin(3) + t88;
t98 = -0.2e1 * qJD(4) * t59 - t76 * t24 + t78 * t27;
t15 = (t58 * t65 - t44) * pkin(7) + (t58 * t59 + t53) * pkin(4) + t98;
t103 = 0.2e1 * qJD(4) * t58 + t78 * t24 + t76 * t27;
t16 = -t57 * pkin(4) + t43 * pkin(7) - t65 * t42 + t103;
t29 = -t32 * mrSges(6,1) + t33 * mrSges(6,2);
t51 = qJDD(5) + t53;
t13 = m(6) * (t83 * t15 - t80 * t16) - t23 * mrSges(6,3) + t51 * mrSges(6,1) - t33 * t29 + t63 * t30;
t14 = m(6) * (t80 * t15 + t83 * t16) + t22 * mrSges(6,3) - t51 * mrSges(6,2) + t32 * t29 - t63 * t31;
t34 = -t58 * mrSges(5,1) + t59 * mrSges(5,2);
t10 = m(5) * t98 + t53 * mrSges(5,1) - t44 * mrSges(5,3) + t83 * t13 + t80 * t14 - t59 * t34 + t65 * t40;
t11 = m(5) * t103 - t53 * mrSges(5,2) + t43 * mrSges(5,3) - t80 * t13 + t83 * t14 + t58 * t34 - t65 * t41;
t61 = qJD(3) * mrSges(4,1) - t66 * mrSges(4,3);
t7 = m(4) * t109 - qJDD(3) * mrSges(4,2) - t53 * mrSges(4,3) - qJD(3) * t61 - t76 * t10 + t78 * t11 - t65 * t49;
t94 = -t79 * mrSges(3,1) + t77 * mrSges(3,2);
t93 = qJDD(1) * mrSges(3,3) + t86 * t94;
t4 = m(3) * t100 + t81 * t7 + t111 * t12 + (-m(3) * t67 - t93) * t77;
t5 = m(3) * t99 + t111 * t7 - t81 * t12 + t93 * t79;
t113 = t79 * t4 + t77 * t5;
t90 = m(4) * t88 + t53 * mrSges(4,1) + t54 * mrSges(4,2) + t78 * t10 + t76 * t11 + t65 * t60 + t66 * t61;
t89 = m(3) * (-qJDD(1) * pkin(1) - t86 * qJ(2) + t96) + t90;
t6 = m(2) * t101 + (-mrSges(2,2) + t114) * t86 + (mrSges(2,1) - t94) * qJDD(1) - t89;
t1 = m(2) * t97 - t86 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t77 * t4 + t79 * t5;
t2 = [-m(1) * g(1) + t84 * t1 - t82 * t6, t1, t5, t7, t11, t14; -m(1) * g(2) + t82 * t1 + t84 * t6, t6, t4, t12, t10, t13; (-m(1) - m(2)) * g(3) + t113, -m(2) * g(3) + t113, t94 * qJDD(1) - t86 * t114 + t89, t90, t87, -t91;];
f_new = t2;
