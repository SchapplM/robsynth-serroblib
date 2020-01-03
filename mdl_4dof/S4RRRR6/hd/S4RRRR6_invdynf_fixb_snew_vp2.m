% Calculate vector of cutting forces with Newton-Euler
% S4RRRR6
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% f_new [3x5]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4RRRR6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(8,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR6_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR6_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR6_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR6_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR6_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:29:25
% EndTime: 2019-12-31 17:29:28
% DurationCPUTime: 1.08s
% Computational Cost: add. (11047->132), mult. (23960->190), div. (0->0), fcn. (17772->10), ass. (0->76)
t58 = sin(pkin(4));
t62 = sin(qJ(2));
t66 = cos(qJ(2));
t78 = qJD(1) * qJD(2);
t50 = (-qJDD(1) * t66 + t62 * t78) * t58;
t88 = pkin(6) * t58;
t59 = cos(pkin(4));
t87 = t59 * g(3);
t68 = qJD(1) ^ 2;
t63 = sin(qJ(1));
t67 = cos(qJ(1));
t75 = t63 * g(1) - t67 * g(2);
t45 = qJDD(1) * pkin(1) + t68 * t88 + t75;
t86 = t45 * t59;
t85 = t58 * t62;
t84 = t58 * t66;
t81 = qJD(1) * t58;
t48 = (-pkin(2) * t66 - pkin(7) * t62) * t81;
t55 = t59 * qJD(1) + qJD(2);
t53 = t55 ^ 2;
t54 = t59 * qJDD(1) + qJDD(2);
t80 = qJD(1) * t66;
t74 = -t67 * g(1) - t63 * g(2);
t46 = -t68 * pkin(1) + qJDD(1) * t88 + t74;
t82 = t66 * t46 + t62 * t86;
t21 = -t53 * pkin(2) + t54 * pkin(7) + (-g(3) * t62 + t48 * t80) * t58 + t82;
t49 = (qJDD(1) * t62 + t66 * t78) * t58;
t22 = t50 * pkin(2) - t49 * pkin(7) - t87 + (-t45 + (pkin(2) * t62 - pkin(7) * t66) * t55 * qJD(1)) * t58;
t61 = sin(qJ(3));
t65 = cos(qJ(3));
t83 = t65 * t21 + t61 * t22;
t77 = t62 * t81;
t38 = t65 * t55 - t61 * t77;
t28 = t38 * qJD(3) + t65 * t49 + t61 * t54;
t39 = t61 * t55 + t65 * t77;
t29 = -t38 * mrSges(4,1) + t39 * mrSges(4,2);
t76 = t58 * t80;
t52 = qJD(3) - t76;
t33 = -t52 * mrSges(4,2) + t38 * mrSges(4,3);
t42 = qJDD(3) + t50;
t60 = sin(qJ(4));
t64 = cos(qJ(4));
t32 = t64 * t39 + t60 * t52;
t16 = -t32 * qJD(4) - t60 * t28 + t64 * t42;
t31 = -t60 * t39 + t64 * t52;
t17 = t31 * qJD(4) + t64 * t28 + t60 * t42;
t37 = qJD(4) - t38;
t24 = -t37 * mrSges(5,2) + t31 * mrSges(5,3);
t25 = t37 * mrSges(5,1) - t32 * mrSges(5,3);
t30 = -t38 * pkin(3) - t39 * pkin(8);
t51 = t52 ^ 2;
t73 = -t61 * t21 + t65 * t22;
t70 = m(5) * (-t42 * pkin(3) - t51 * pkin(8) + t39 * t30 - t73) - t16 * mrSges(5,1) + t17 * mrSges(5,2) - t31 * t24 + t32 * t25;
t10 = m(4) * t73 + t42 * mrSges(4,1) - t28 * mrSges(4,3) - t39 * t29 + t52 * t33 - t70;
t43 = t55 * mrSges(3,1) - mrSges(3,3) * t77;
t47 = (-mrSges(3,1) * t66 + mrSges(3,2) * t62) * t81;
t14 = -t51 * pkin(3) + t42 * pkin(8) + t38 * t30 + t83;
t72 = -g(3) * t84 - t62 * t46 + t66 * t86;
t20 = -t54 * pkin(2) - t53 * pkin(7) + t48 * t77 - t72;
t27 = -t39 * qJD(3) - t61 * t49 + t65 * t54;
t15 = (-t38 * t52 - t28) * pkin(8) + (t39 * t52 - t27) * pkin(3) + t20;
t23 = -t31 * mrSges(5,1) + t32 * mrSges(5,2);
t26 = qJDD(4) - t27;
t11 = m(5) * (-t60 * t14 + t64 * t15) - t17 * mrSges(5,3) + t26 * mrSges(5,1) - t32 * t23 + t37 * t24;
t12 = m(5) * (t64 * t14 + t60 * t15) + t16 * mrSges(5,3) - t26 * mrSges(5,2) + t31 * t23 - t37 * t25;
t34 = t52 * mrSges(4,1) - t39 * mrSges(4,3);
t9 = m(4) * t83 - t42 * mrSges(4,2) + t27 * mrSges(4,3) - t60 * t11 + t64 * t12 + t38 * t29 - t52 * t34;
t4 = m(3) * (-g(3) * t85 + t82) - t50 * mrSges(3,3) - t54 * mrSges(3,2) + t47 * t76 - t55 * t43 + t65 * t9 - t61 * t10;
t44 = -t55 * mrSges(3,2) + mrSges(3,3) * t76;
t6 = m(3) * (-t58 * t45 - t87) + t49 * mrSges(3,2) + t50 * mrSges(3,1) + t61 * t9 + t65 * t10 + (t43 * t62 - t44 * t66) * t81;
t69 = m(4) * t20 - t27 * mrSges(4,1) + t28 * mrSges(4,2) + t64 * t11 + t60 * t12 - t38 * t33 + t39 * t34;
t8 = m(3) * t72 + t54 * mrSges(3,1) - t49 * mrSges(3,3) + t55 * t44 - t47 * t77 - t69;
t79 = t4 * t85 + t59 * t6 + t8 * t84;
t2 = m(2) * t74 - t68 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t66 * t4 - t62 * t8;
t1 = m(2) * t75 + qJDD(1) * mrSges(2,1) - t68 * mrSges(2,2) - t58 * t6 + (t62 * t4 + t66 * t8) * t59;
t3 = [-m(1) * g(1) - t63 * t1 + t67 * t2, t2, t4, t9, t12; -m(1) * g(2) + t67 * t1 + t63 * t2, t1, t8, t10, t11; (-m(1) - m(2)) * g(3) + t79, -m(2) * g(3) + t79, t6, t69, t70;];
f_new = t3;
