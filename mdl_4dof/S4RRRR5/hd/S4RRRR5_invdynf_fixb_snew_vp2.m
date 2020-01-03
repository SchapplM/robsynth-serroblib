% Calculate vector of cutting forces with Newton-Euler
% S4RRRR5
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
% Datum: 2019-12-31 17:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4RRRR5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR5_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR5_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR5_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR5_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR5_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR5_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR5_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:27:24
% EndTime: 2019-12-31 17:27:26
% DurationCPUTime: 0.68s
% Computational Cost: add. (6467->123), mult. (12902->164), div. (0->0), fcn. (8203->8), ass. (0->67)
t66 = qJD(1) ^ 2;
t60 = sin(qJ(1));
t64 = cos(qJ(1));
t75 = t60 * g(1) - t64 * g(2);
t39 = -qJDD(1) * pkin(1) - t66 * pkin(5) - t75;
t59 = sin(qJ(2));
t63 = cos(qJ(2));
t76 = qJD(1) * qJD(2);
t73 = t63 * t76;
t48 = t59 * qJDD(1) + t73;
t54 = t59 * t76;
t49 = t63 * qJDD(1) - t54;
t58 = sin(qJ(3));
t62 = cos(qJ(3));
t78 = qJD(1) * t59;
t44 = t62 * qJD(2) - t58 * t78;
t29 = t44 * qJD(3) + t58 * qJDD(2) + t62 * t48;
t43 = qJDD(3) - t49;
t45 = t58 * qJD(2) + t62 * t78;
t77 = t63 * qJD(1);
t53 = qJD(3) - t77;
t22 = (-t48 - t73) * pkin(6) + (-t49 + t54) * pkin(2) + t39;
t47 = (-pkin(2) * t63 - pkin(6) * t59) * qJD(1);
t65 = qJD(2) ^ 2;
t71 = -t64 * g(1) - t60 * g(2);
t40 = -t66 * pkin(1) + qJDD(1) * pkin(5) + t71;
t74 = -t59 * g(3) + t63 * t40;
t25 = -t65 * pkin(2) + qJDD(2) * pkin(6) + t47 * t77 + t74;
t72 = t62 * t22 - t58 * t25;
t11 = (t44 * t53 - t29) * pkin(7) + (t44 * t45 + t43) * pkin(3) + t72;
t28 = -t45 * qJD(3) + t62 * qJDD(2) - t58 * t48;
t35 = t53 * pkin(3) - t45 * pkin(7);
t42 = t44 ^ 2;
t80 = t58 * t22 + t62 * t25;
t12 = -t42 * pkin(3) + t28 * pkin(7) - t53 * t35 + t80;
t57 = sin(qJ(4));
t61 = cos(qJ(4));
t31 = t57 * t44 + t61 * t45;
t16 = -t31 * qJD(4) + t61 * t28 - t57 * t29;
t30 = t61 * t44 - t57 * t45;
t19 = -t30 * mrSges(5,1) + t31 * mrSges(5,2);
t52 = qJD(4) + t53;
t27 = t52 * mrSges(5,1) - t31 * mrSges(5,3);
t41 = qJDD(4) + t43;
t10 = m(5) * (t57 * t11 + t61 * t12) + t16 * mrSges(5,3) - t41 * mrSges(5,2) + t30 * t19 - t52 * t27;
t32 = -t44 * mrSges(4,1) + t45 * mrSges(4,2);
t33 = -t53 * mrSges(4,2) + t44 * mrSges(4,3);
t17 = t30 * qJD(4) + t57 * t28 + t61 * t29;
t26 = -t52 * mrSges(5,2) + t30 * mrSges(5,3);
t9 = m(5) * (t61 * t11 - t57 * t12) - t17 * mrSges(5,3) + t41 * mrSges(5,1) - t31 * t19 + t52 * t26;
t5 = m(4) * t72 + t43 * mrSges(4,1) - t29 * mrSges(4,3) + t57 * t10 - t45 * t32 + t53 * t33 + t61 * t9;
t50 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t78;
t51 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t77;
t34 = t53 * mrSges(4,1) - t45 * mrSges(4,3);
t6 = m(4) * t80 - t43 * mrSges(4,2) + t28 * mrSges(4,3) + t61 * t10 + t44 * t32 - t53 * t34 - t57 * t9;
t82 = m(3) * t39 - t49 * mrSges(3,1) + t48 * mrSges(3,2) + t62 * t5 + t58 * t6 + (t59 * t50 - t63 * t51) * qJD(1);
t46 = (-mrSges(3,1) * t63 + mrSges(3,2) * t59) * qJD(1);
t4 = m(3) * t74 - qJDD(2) * mrSges(3,2) + t49 * mrSges(3,3) - qJD(2) * t50 + t46 * t77 - t58 * t5 + t62 * t6;
t79 = -t63 * g(3) - t59 * t40;
t24 = -qJDD(2) * pkin(2) - t65 * pkin(6) + t47 * t78 - t79;
t69 = t16 * mrSges(5,1) + t30 * t26 - m(5) * (-t28 * pkin(3) - t42 * pkin(7) + t45 * t35 + t24) - t17 * mrSges(5,2) - t31 * t27;
t67 = m(4) * t24 - t28 * mrSges(4,1) + t29 * mrSges(4,2) - t44 * t33 + t45 * t34 - t69;
t8 = m(3) * t79 + qJDD(2) * mrSges(3,1) - t48 * mrSges(3,3) + qJD(2) * t51 - t46 * t78 - t67;
t81 = t59 * t4 + t63 * t8;
t2 = m(2) * t75 + qJDD(1) * mrSges(2,1) - t66 * mrSges(2,2) - t82;
t1 = m(2) * t71 - t66 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t63 * t4 - t59 * t8;
t3 = [-m(1) * g(1) + t64 * t1 - t60 * t2, t1, t4, t6, t10; -m(1) * g(2) + t60 * t1 + t64 * t2, t2, t8, t5, t9; (-m(1) - m(2)) * g(3) + t81, -m(2) * g(3) + t81, t82, t67, -t69;];
f_new = t3;
