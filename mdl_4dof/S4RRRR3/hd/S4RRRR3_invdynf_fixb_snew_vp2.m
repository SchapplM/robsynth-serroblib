% Calculate vector of cutting forces with Newton-Euler
% S4RRRR3
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 17:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4RRRR3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR3_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR3_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR3_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR3_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR3_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR3_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR3_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:24:18
% EndTime: 2019-12-31 17:24:20
% DurationCPUTime: 0.77s
% Computational Cost: add. (7517->126), mult. (16244->170), div. (0->0), fcn. (10653->8), ass. (0->66)
t60 = sin(qJ(2));
t64 = cos(qJ(2));
t77 = qJD(1) * qJD(2);
t47 = t60 * qJDD(1) + t64 * t77;
t48 = t64 * qJDD(1) - t60 * t77;
t79 = qJD(1) * t60;
t49 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t79;
t78 = qJD(1) * t64;
t50 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t78;
t66 = qJD(1) ^ 2;
t59 = sin(qJ(3));
t63 = cos(qJ(3));
t42 = (t59 * t64 + t60 * t63) * qJD(1);
t27 = -t42 * qJD(3) - t59 * t47 + t63 * t48;
t41 = (-t59 * t60 + t63 * t64) * qJD(1);
t28 = t41 * qJD(3) + t63 * t47 + t59 * t48;
t56 = qJD(2) + qJD(3);
t36 = -t56 * mrSges(4,2) + t41 * mrSges(4,3);
t37 = t56 * mrSges(4,1) - t42 * mrSges(4,3);
t51 = qJD(2) * pkin(2) - pkin(6) * t79;
t57 = t64 ^ 2;
t61 = sin(qJ(1));
t65 = cos(qJ(1));
t76 = t61 * g(1) - t65 * g(2);
t71 = -qJDD(1) * pkin(1) - t76;
t69 = -t48 * pkin(2) + t51 * t79 + (-pkin(6) * t57 - pkin(5)) * t66 + t71;
t58 = sin(qJ(4));
t62 = cos(qJ(4));
t34 = t58 * t41 + t62 * t42;
t16 = -t34 * qJD(4) + t62 * t27 - t58 * t28;
t33 = t62 * t41 - t58 * t42;
t17 = t33 * qJD(4) + t58 * t27 + t62 * t28;
t53 = qJD(4) + t56;
t30 = -t53 * mrSges(5,2) + t33 * mrSges(5,3);
t31 = t53 * mrSges(5,1) - t34 * mrSges(5,3);
t38 = t56 * pkin(3) - t42 * pkin(7);
t40 = t41 ^ 2;
t70 = t16 * mrSges(5,1) + t33 * t30 - m(5) * (-t27 * pkin(3) - t40 * pkin(7) + t42 * t38 + t69) - t17 * mrSges(5,2) - t34 * t31;
t68 = -m(4) * t69 + t27 * mrSges(4,1) - t28 * mrSges(4,2) + t41 * t36 - t42 * t37 + t70;
t84 = (t60 * t49 - t64 * t50) * qJD(1) + m(3) * (-t66 * pkin(5) + t71) - t48 * mrSges(3,1) + t47 * mrSges(3,2) - t68;
t46 = (-mrSges(3,1) * t64 + mrSges(3,2) * t60) * qJD(1);
t55 = qJDD(2) + qJDD(3);
t73 = -t65 * g(1) - t61 * g(2);
t44 = -t66 * pkin(1) + qJDD(1) * pkin(5) + t73;
t81 = t60 * t44;
t82 = pkin(2) * t66;
t23 = qJDD(2) * pkin(2) - t47 * pkin(6) - t81 + (pkin(6) * t77 + t60 * t82 - g(3)) * t64;
t75 = -t60 * g(3) + t64 * t44;
t24 = t48 * pkin(6) - qJD(2) * t51 - t57 * t82 + t75;
t74 = t63 * t23 - t59 * t24;
t11 = (t41 * t56 - t28) * pkin(7) + (t41 * t42 + t55) * pkin(3) + t74;
t80 = t59 * t23 + t63 * t24;
t12 = -t40 * pkin(3) + t27 * pkin(7) - t56 * t38 + t80;
t19 = -t33 * mrSges(5,1) + t34 * mrSges(5,2);
t52 = qJDD(4) + t55;
t10 = m(5) * (t58 * t11 + t62 * t12) + t16 * mrSges(5,3) - t52 * mrSges(5,2) + t33 * t19 - t53 * t31;
t35 = -t41 * mrSges(4,1) + t42 * mrSges(4,2);
t9 = m(5) * (t62 * t11 - t58 * t12) - t17 * mrSges(5,3) + t52 * mrSges(5,1) - t34 * t19 + t53 * t30;
t6 = m(4) * t74 + t55 * mrSges(4,1) - t28 * mrSges(4,3) + t58 * t10 - t42 * t35 + t56 * t36 + t62 * t9;
t7 = m(4) * t80 - t55 * mrSges(4,2) + t27 * mrSges(4,3) + t62 * t10 + t41 * t35 - t56 * t37 - t58 * t9;
t4 = m(3) * (-t64 * g(3) - t81) - t47 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t46 * t79 + qJD(2) * t50 + t59 * t7 + t63 * t6;
t5 = m(3) * t75 - qJDD(2) * mrSges(3,2) + t48 * mrSges(3,3) - qJD(2) * t49 + t46 * t78 - t59 * t6 + t63 * t7;
t83 = t64 * t4 + t60 * t5;
t8 = m(2) * t76 + qJDD(1) * mrSges(2,1) - t66 * mrSges(2,2) - t84;
t1 = m(2) * t73 - t66 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t60 * t4 + t64 * t5;
t2 = [-m(1) * g(1) + t65 * t1 - t61 * t8, t1, t5, t7, t10; -m(1) * g(2) + t61 * t1 + t65 * t8, t8, t4, t6, t9; (-m(1) - m(2)) * g(3) + t83, -m(2) * g(3) + t83, t84, -t68, -t70;];
f_new = t2;
