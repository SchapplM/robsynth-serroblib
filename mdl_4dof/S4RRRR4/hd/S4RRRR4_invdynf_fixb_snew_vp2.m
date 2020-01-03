% Calculate vector of cutting forces with Newton-Euler
% S4RRRR4
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
% Datum: 2019-12-31 17:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4RRRR4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR4_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR4_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR4_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR4_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR4_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR4_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR4_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:25:42
% EndTime: 2019-12-31 17:25:44
% DurationCPUTime: 0.67s
% Computational Cost: add. (5981->125), mult. (12234->170), div. (0->0), fcn. (7881->8), ass. (0->66)
t57 = sin(qJ(3));
t58 = sin(qJ(2));
t61 = cos(qJ(3));
t62 = cos(qJ(2));
t40 = (t57 * t58 - t61 * t62) * qJD(1);
t76 = qJD(1) * qJD(2);
t46 = t58 * qJDD(1) + t62 * t76;
t47 = t62 * qJDD(1) - t58 * t76;
t78 = qJD(1) * t58;
t48 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t78;
t77 = qJD(1) * t62;
t49 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t77;
t64 = qJD(1) ^ 2;
t41 = (t57 * t62 + t58 * t61) * qJD(1);
t26 = -t41 * qJD(3) - t57 * t46 + t61 * t47;
t27 = -t40 * qJD(3) + t61 * t46 + t57 * t47;
t54 = qJD(2) + qJD(3);
t50 = qJD(2) * pkin(2) - pkin(6) * t78;
t55 = t62 ^ 2;
t59 = sin(qJ(1));
t63 = cos(qJ(1));
t75 = t59 * g(1) - t63 * g(2);
t69 = -qJDD(1) * pkin(1) - t75;
t66 = -t47 * pkin(2) + t50 * t78 + (-pkin(6) * t55 - pkin(5)) * t64 + t69;
t13 = (t40 * t54 - t27) * pkin(7) + (t41 * t54 - t26) * pkin(3) + t66;
t33 = t40 * pkin(3) - t41 * pkin(7);
t52 = t54 ^ 2;
t53 = qJDD(2) + qJDD(3);
t73 = -t63 * g(1) - t59 * g(2);
t43 = -t64 * pkin(1) + qJDD(1) * pkin(5) + t73;
t80 = t58 * t43;
t81 = pkin(2) * t64;
t21 = qJDD(2) * pkin(2) - t46 * pkin(6) - t80 + (pkin(6) * t76 + t58 * t81 - g(3)) * t62;
t74 = -t58 * g(3) + t62 * t43;
t22 = t47 * pkin(6) - qJD(2) * t50 - t55 * t81 + t74;
t79 = t57 * t21 + t61 * t22;
t15 = -t52 * pkin(3) + t53 * pkin(7) - t40 * t33 + t79;
t56 = sin(qJ(4));
t60 = cos(qJ(4));
t34 = -t56 * t41 + t60 * t54;
t17 = t34 * qJD(4) + t60 * t27 + t56 * t53;
t35 = t60 * t41 + t56 * t54;
t20 = -t34 * mrSges(5,1) + t35 * mrSges(5,2);
t24 = qJDD(4) - t26;
t39 = qJD(4) + t40;
t29 = -t39 * mrSges(5,2) + t34 * mrSges(5,3);
t11 = m(5) * (t60 * t13 - t56 * t15) - t17 * mrSges(5,3) + t24 * mrSges(5,1) - t35 * t20 + t39 * t29;
t16 = -t35 * qJD(4) - t56 * t27 + t60 * t53;
t30 = t39 * mrSges(5,1) - t35 * mrSges(5,3);
t12 = m(5) * (t56 * t13 + t60 * t15) + t16 * mrSges(5,3) - t24 * mrSges(5,2) + t34 * t20 - t39 * t30;
t36 = -t54 * mrSges(4,2) - t40 * mrSges(4,3);
t37 = t54 * mrSges(4,1) - t41 * mrSges(4,3);
t68 = -m(4) * t66 + t26 * mrSges(4,1) - t27 * mrSges(4,2) - t60 * t11 - t56 * t12 - t40 * t36 - t41 * t37;
t83 = (t58 * t48 - t62 * t49) * qJD(1) + m(3) * (-t64 * pkin(5) + t69) - t47 * mrSges(3,1) + t46 * mrSges(3,2) - t68;
t45 = (-mrSges(3,1) * t62 + mrSges(3,2) * t58) * qJD(1);
t32 = t40 * mrSges(4,1) + t41 * mrSges(4,2);
t7 = m(4) * t79 - t53 * mrSges(4,2) + t26 * mrSges(4,3) - t56 * t11 + t60 * t12 - t40 * t32 - t54 * t37;
t72 = t61 * t21 - t57 * t22;
t67 = m(5) * (-t53 * pkin(3) - t52 * pkin(7) + t41 * t33 - t72) - t16 * mrSges(5,1) + t17 * mrSges(5,2) - t34 * t29 + t35 * t30;
t8 = m(4) * t72 + t53 * mrSges(4,1) - t27 * mrSges(4,3) - t41 * t32 + t54 * t36 - t67;
t4 = m(3) * (-t62 * g(3) - t80) - t46 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t45 * t78 + qJD(2) * t49 + t57 * t7 + t61 * t8;
t5 = m(3) * t74 - qJDD(2) * mrSges(3,2) + t47 * mrSges(3,3) - qJD(2) * t48 + t45 * t77 - t57 * t8 + t61 * t7;
t82 = t62 * t4 + t58 * t5;
t6 = m(2) * t75 + qJDD(1) * mrSges(2,1) - t64 * mrSges(2,2) - t83;
t1 = m(2) * t73 - t64 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t58 * t4 + t62 * t5;
t2 = [-m(1) * g(1) + t63 * t1 - t59 * t6, t1, t5, t7, t12; -m(1) * g(2) + t59 * t1 + t63 * t6, t6, t4, t8, t11; (-m(1) - m(2)) * g(3) + t82, -m(2) * g(3) + t82, t83, -t68, t67;];
f_new = t2;
