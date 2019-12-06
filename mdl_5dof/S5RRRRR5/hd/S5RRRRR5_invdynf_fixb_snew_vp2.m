% Calculate vector of cutting forces with Newton-Euler
% S5RRRRR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-05 18:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRRR5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR5_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR5_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR5_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR5_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:58:14
% EndTime: 2019-12-05 18:58:16
% DurationCPUTime: 0.80s
% Computational Cost: add. (14593->116), mult. (15052->154), div. (0->0), fcn. (8503->10), ass. (0->69)
t53 = qJDD(1) + qJDD(2);
t46 = qJDD(3) + t53;
t58 = sin(qJ(4));
t63 = cos(qJ(4));
t55 = qJD(1) + qJD(2);
t47 = qJD(3) + t55;
t77 = qJD(4) * t47;
t35 = t58 * t46 + t63 * t77;
t36 = t63 * t46 - t58 * t77;
t83 = t47 * t58;
t40 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t83;
t82 = t47 * t63;
t41 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t82;
t45 = t47 ^ 2;
t57 = sin(qJ(5));
t62 = cos(qJ(5));
t33 = (t57 * t63 + t58 * t62) * t47;
t21 = -t33 * qJD(5) - t57 * t35 + t62 * t36;
t32 = (-t57 * t58 + t62 * t63) * t47;
t22 = t32 * qJD(5) + t62 * t35 + t57 * t36;
t54 = qJD(4) + qJD(5);
t30 = -t54 * mrSges(6,2) + t32 * mrSges(6,3);
t31 = t54 * mrSges(6,1) - t33 * mrSges(6,3);
t42 = qJD(4) * pkin(4) - pkin(9) * t83;
t56 = t63 ^ 2;
t61 = sin(qJ(1));
t66 = cos(qJ(1));
t78 = t66 * g(2) + t61 * g(3);
t43 = qJDD(1) * pkin(1) + t78;
t67 = qJD(1) ^ 2;
t74 = t61 * g(2) - t66 * g(3);
t44 = -t67 * pkin(1) + t74;
t60 = sin(qJ(2));
t65 = cos(qJ(2));
t72 = t65 * t43 - t60 * t44;
t28 = t53 * pkin(2) + t72;
t51 = t55 ^ 2;
t79 = t60 * t43 + t65 * t44;
t29 = -t51 * pkin(2) + t79;
t59 = sin(qJ(3));
t64 = cos(qJ(3));
t73 = t64 * t28 - t59 * t29;
t70 = -t46 * pkin(3) - t73;
t69 = t21 * mrSges(6,1) + t32 * t30 - m(6) * (t42 * t83 - t36 * pkin(4) + (-pkin(9) * t56 - pkin(8)) * t45 + t70) - t22 * mrSges(6,2) - t33 * t31;
t87 = (t58 * t40 - t63 * t41) * t47 + m(5) * (-t45 * pkin(8) + t70) - t36 * mrSges(5,1) + t35 * mrSges(5,2) - t69;
t86 = -m(3) - m(4);
t80 = t59 * t28 + t64 * t29;
t19 = -t45 * pkin(3) + t46 * pkin(8) + t80;
t81 = t58 * t19;
t84 = pkin(4) * t45;
t13 = qJDD(4) * pkin(4) - t35 * pkin(9) - t81 + (pkin(9) * t77 + t58 * t84 - g(1)) * t63;
t75 = -t58 * g(1) + t63 * t19;
t14 = t36 * pkin(9) - qJD(4) * t42 - t56 * t84 + t75;
t24 = -t32 * mrSges(6,1) + t33 * mrSges(6,2);
t52 = qJDD(4) + qJDD(5);
t11 = m(6) * (t62 * t13 - t57 * t14) - t22 * mrSges(6,3) + t52 * mrSges(6,1) - t33 * t24 + t54 * t30;
t12 = m(6) * (t57 * t13 + t62 * t14) + t21 * mrSges(6,3) - t52 * mrSges(6,2) + t32 * t24 - t54 * t31;
t34 = (-mrSges(5,1) * t63 + mrSges(5,2) * t58) * t47;
t8 = m(5) * (-t63 * g(1) - t81) - t35 * mrSges(5,3) + qJDD(4) * mrSges(5,1) - t34 * t83 + qJD(4) * t41 + t57 * t12 + t62 * t11;
t9 = m(5) * t75 - qJDD(4) * mrSges(5,2) + t36 * mrSges(5,3) - qJD(4) * t40 - t57 * t11 + t62 * t12 + t34 * t82;
t85 = t58 * t9 + t63 * t8;
t76 = -m(2) + t86;
t10 = m(4) * t73 + t46 * mrSges(4,1) - t45 * mrSges(4,2) - t87;
t5 = m(4) * t80 - t45 * mrSges(4,1) - t46 * mrSges(4,2) - t58 * t8 + t63 * t9;
t4 = m(3) * t79 - t51 * mrSges(3,1) - t53 * mrSges(3,2) - t59 * t10 + t64 * t5;
t3 = m(3) * t72 + t53 * mrSges(3,1) - t51 * mrSges(3,2) + t64 * t10 + t59 * t5;
t2 = m(2) * t78 + qJDD(1) * mrSges(2,1) - t67 * mrSges(2,2) + t65 * t3 + t60 * t4;
t1 = m(2) * t74 - t67 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t60 * t3 + t65 * t4;
t6 = [(-m(1) + t76) * g(1) + t85, t1, t4, t5, t9, t12; -m(1) * g(2) - t61 * t1 - t66 * t2, t2, t3, t10, t8, t11; -m(1) * g(3) + t66 * t1 - t61 * t2, t76 * g(1) + t85, t86 * g(1) + t85, -m(4) * g(1) + t85, t87, -t69;];
f_new = t6;
