% Calculate vector of cutting forces with Newton-Euler
% S5RPRRP2
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% m [6x1]
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
% Datum: 2022-01-23 09:28
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRRP2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP2_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP2_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP2_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP2_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:27:54
% EndTime: 2022-01-23 09:27:55
% DurationCPUTime: 0.52s
% Computational Cost: add. (5154->109), mult. (7066->137), div. (0->0), fcn. (3444->8), ass. (0->58)
t48 = qJDD(1) + qJDD(3);
t55 = sin(qJ(4));
t58 = cos(qJ(4));
t49 = qJD(1) + qJD(3);
t75 = qJD(4) * t49;
t69 = t58 * t75;
t29 = t55 * t48 + t69;
t30 = t58 * t48 - t55 * t75;
t83 = t49 * t55;
t39 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t83;
t82 = t49 * t58;
t40 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t82;
t41 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t82;
t47 = t49 ^ 2;
t57 = sin(qJ(1));
t60 = cos(qJ(1));
t68 = t57 * g(1) - t60 * g(2);
t35 = qJDD(1) * pkin(1) + t68;
t61 = qJD(1) ^ 2;
t64 = -t60 * g(1) - t57 * g(2);
t36 = -t61 * pkin(1) + t64;
t53 = sin(pkin(8));
t54 = cos(pkin(8));
t65 = t54 * t35 - t53 * t36;
t21 = qJDD(1) * pkin(2) + t65;
t77 = t53 * t35 + t54 * t36;
t22 = -t61 * pkin(2) + t77;
t56 = sin(qJ(3));
t59 = cos(qJ(3));
t66 = t59 * t21 - t56 * t22;
t63 = -t48 * pkin(3) - t66;
t37 = qJD(4) * pkin(4) - qJ(5) * t83;
t38 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t83;
t51 = t58 ^ 2;
t70 = m(6) * (t37 * t83 - t30 * pkin(4) + qJDD(5) + (-qJ(5) * t51 - pkin(7)) * t47 + t63) + t38 * t83 + t29 * mrSges(6,2);
t87 = -(-t55 * t39 + (t40 + t41) * t58) * t49 - (mrSges(5,1) + mrSges(6,1)) * t30 + m(5) * (-t47 * pkin(7) + t63) + t29 * mrSges(5,2) + t70;
t84 = pkin(4) * t47;
t78 = t56 * t21 + t59 * t22;
t17 = -t47 * pkin(3) + t48 * pkin(7) + t78;
t52 = -g(3) + qJDD(2);
t79 = t58 * t17 + t55 * t52;
t74 = qJD(5) * t49;
t28 = (-mrSges(5,1) * t58 + mrSges(5,2) * t55) * t49;
t27 = (-mrSges(6,1) * t58 + mrSges(6,2) * t55) * t49;
t71 = m(6) * (t30 * qJ(5) - qJD(4) * t37 - t51 * t84 + 0.2e1 * t58 * t74 + t79) + t27 * t82 + t30 * mrSges(6,3);
t10 = m(5) * t79 + t30 * mrSges(5,3) + t28 * t82 + (-mrSges(5,2) - mrSges(6,2)) * qJDD(4) + (-t39 - t38) * qJD(4) + t71;
t43 = t58 * t52;
t72 = m(6) * (qJDD(4) * pkin(4) + t43 + (-t29 + t69) * qJ(5) + (t58 * t84 - t17 - 0.2e1 * t74) * t55) + qJD(4) * t40 + qJDD(4) * mrSges(6,1);
t9 = m(5) * (-t55 * t17 + t43) + qJDD(4) * mrSges(5,1) + qJD(4) * t41 + (-t27 - t28) * t83 + (-mrSges(5,3) - mrSges(6,3)) * t29 + t72;
t73 = m(4) * t52 + t55 * t10 + t58 * t9;
t67 = m(3) * t52 + t73;
t6 = m(4) * t66 + t48 * mrSges(4,1) - t47 * mrSges(4,2) - t87;
t5 = m(4) * t78 - t47 * mrSges(4,1) - t48 * mrSges(4,2) + t58 * t10 - t55 * t9;
t4 = m(3) * t77 - t61 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t59 * t5 - t56 * t6;
t3 = m(3) * t65 + qJDD(1) * mrSges(3,1) - t61 * mrSges(3,2) + t56 * t5 + t59 * t6;
t2 = m(2) * t64 - t61 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t53 * t3 + t54 * t4;
t1 = m(2) * t68 + qJDD(1) * mrSges(2,1) - t61 * mrSges(2,2) + t54 * t3 + t53 * t4;
t7 = [-m(1) * g(1) - t57 * t1 + t60 * t2, t2, t4, t5, t10, -qJDD(4) * mrSges(6,2) - qJD(4) * t38 + t71; -m(1) * g(2) + t60 * t1 + t57 * t2, t1, t3, t6, t9, -t29 * mrSges(6,3) - t27 * t83 + t72; (-m(1) - m(2)) * g(3) + t67, -m(2) * g(3) + t67, t67, t73, t87, -t30 * mrSges(6,1) - t40 * t82 + t70;];
f_new = t7;
