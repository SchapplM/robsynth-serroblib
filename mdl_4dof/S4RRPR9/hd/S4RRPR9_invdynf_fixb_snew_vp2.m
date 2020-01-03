% Calculate vector of cutting forces with Newton-Euler
% S4RRPR9
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4RRPR9_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR9_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR9_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR9_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR9_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR9_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR9_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR9_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR9_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:09:23
% EndTime: 2019-12-31 17:09:24
% DurationCPUTime: 0.64s
% Computational Cost: add. (5764->122), mult. (12432->166), div. (0->0), fcn. (7765->8), ass. (0->65)
t66 = qJD(1) ^ 2;
t61 = sin(qJ(1));
t64 = cos(qJ(1));
t75 = t61 * g(1) - t64 * g(2);
t41 = -qJDD(1) * pkin(1) - t66 * pkin(5) - t75;
t60 = sin(qJ(2));
t63 = cos(qJ(2));
t77 = qJD(1) * qJD(2);
t73 = t63 * t77;
t49 = t60 * qJDD(1) + t73;
t57 = sin(pkin(7));
t58 = cos(pkin(7));
t34 = t57 * qJDD(2) + t58 * t49;
t79 = qJD(1) * t60;
t44 = t58 * qJD(2) - t57 * t79;
t45 = t57 * qJD(2) + t58 * t79;
t54 = t60 * t77;
t50 = t63 * qJDD(1) - t54;
t22 = (-t49 - t73) * qJ(3) + (-t50 + t54) * pkin(2) + t41;
t47 = (-pkin(2) * t63 - qJ(3) * t60) * qJD(1);
t65 = qJD(2) ^ 2;
t71 = -t64 * g(1) - t61 * g(2);
t42 = -t66 * pkin(1) + qJDD(1) * pkin(5) + t71;
t74 = -t60 * g(3) + t63 * t42;
t78 = t63 * qJD(1);
t25 = -t65 * pkin(2) + qJDD(2) * qJ(3) + t47 * t78 + t74;
t72 = -0.2e1 * qJD(3) * t45 + t58 * t22 - t57 * t25;
t11 = (-t44 * t78 - t34) * pkin(6) + (t44 * t45 - t50) * pkin(3) + t72;
t33 = t58 * qJDD(2) - t57 * t49;
t35 = -pkin(3) * t78 - t45 * pkin(6);
t43 = t44 ^ 2;
t76 = 0.2e1 * qJD(3) * t44 + t57 * t22 + t58 * t25;
t12 = -t43 * pkin(3) + t33 * pkin(6) + t35 * t78 + t76;
t59 = sin(qJ(4));
t62 = cos(qJ(4));
t29 = t59 * t44 + t62 * t45;
t16 = -t29 * qJD(4) + t62 * t33 - t59 * t34;
t28 = t62 * t44 - t59 * t45;
t19 = -t28 * mrSges(5,1) + t29 * mrSges(5,2);
t53 = qJD(4) - t78;
t27 = t53 * mrSges(5,1) - t29 * mrSges(5,3);
t46 = qJDD(4) - t50;
t10 = m(5) * (t59 * t11 + t62 * t12) + t16 * mrSges(5,3) - t46 * mrSges(5,2) + t28 * t19 - t53 * t27;
t30 = -t44 * mrSges(4,1) + t45 * mrSges(4,2);
t31 = mrSges(4,2) * t78 + t44 * mrSges(4,3);
t17 = t28 * qJD(4) + t59 * t33 + t62 * t34;
t26 = -t53 * mrSges(5,2) + t28 * mrSges(5,3);
t9 = m(5) * (t62 * t11 - t59 * t12) - t17 * mrSges(5,3) + t46 * mrSges(5,1) - t29 * t19 + t53 * t26;
t5 = m(4) * t72 - t50 * mrSges(4,1) - t34 * mrSges(4,3) + t59 * t10 - t45 * t30 - t31 * t78 + t62 * t9;
t51 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t79;
t52 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t78;
t32 = -mrSges(4,1) * t78 - t45 * mrSges(4,3);
t6 = m(4) * t76 + t50 * mrSges(4,2) + t33 * mrSges(4,3) + t62 * t10 + t44 * t30 + t32 * t78 - t59 * t9;
t82 = m(3) * t41 - t50 * mrSges(3,1) + t49 * mrSges(3,2) + t58 * t5 + t57 * t6 + (t60 * t51 - t63 * t52) * qJD(1);
t48 = (-mrSges(3,1) * t63 + mrSges(3,2) * t60) * qJD(1);
t4 = m(3) * t74 - qJDD(2) * mrSges(3,2) + t50 * mrSges(3,3) - qJD(2) * t51 + t48 * t78 - t57 * t5 + t58 * t6;
t80 = -t63 * g(3) - t60 * t42;
t24 = -qJDD(2) * pkin(2) - t65 * qJ(3) + t47 * t79 + qJDD(3) - t80;
t69 = t16 * mrSges(5,1) + t28 * t26 - m(5) * (-t33 * pkin(3) - t43 * pkin(6) + t45 * t35 + t24) - t17 * mrSges(5,2) - t29 * t27;
t67 = m(4) * t24 - t33 * mrSges(4,1) + t34 * mrSges(4,2) - t44 * t31 + t45 * t32 - t69;
t8 = m(3) * t80 + qJDD(2) * mrSges(3,1) - t49 * mrSges(3,3) + qJD(2) * t52 - t48 * t79 - t67;
t81 = t60 * t4 + t63 * t8;
t2 = m(2) * t75 + qJDD(1) * mrSges(2,1) - t66 * mrSges(2,2) - t82;
t1 = m(2) * t71 - t66 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t63 * t4 - t60 * t8;
t3 = [-m(1) * g(1) + t64 * t1 - t61 * t2, t1, t4, t6, t10; -m(1) * g(2) + t61 * t1 + t64 * t2, t2, t8, t5, t9; (-m(1) - m(2)) * g(3) + t81, -m(2) * g(3) + t81, t82, t67, -t69;];
f_new = t3;
