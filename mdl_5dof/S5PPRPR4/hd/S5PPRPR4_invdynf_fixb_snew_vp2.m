% Calculate vector of cutting forces with Newton-Euler
% S5PPRPR4
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PPRPR4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR4_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR4_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR4_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR4_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR4_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR4_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR4_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:32:18
% EndTime: 2019-12-31 17:32:19
% DurationCPUTime: 0.37s
% Computational Cost: add. (3088->82), mult. (6290->112), div. (0->0), fcn. (3995->8), ass. (0->54)
t50 = qJD(3) ^ 2;
t44 = cos(pkin(8));
t40 = t44 ^ 2;
t42 = sin(pkin(8));
t68 = t42 ^ 2 + t40;
t72 = t68 * mrSges(5,3);
t71 = pkin(4) * t50;
t43 = sin(pkin(7));
t45 = cos(pkin(7));
t62 = t43 * g(1) - t45 * g(2);
t30 = qJDD(2) - t62;
t32 = -t45 * g(1) - t43 * g(2);
t47 = sin(qJ(3));
t49 = cos(qJ(3));
t70 = t47 * t30 + t49 * t32;
t41 = g(3) - qJDD(1);
t66 = qJD(3) * qJD(4);
t69 = t44 * t41 - 0.2e1 * t42 * t66;
t67 = pkin(6) * qJDD(3);
t19 = -t50 * pkin(3) + qJDD(3) * qJ(4) + t70;
t65 = t42 * t41 + (t19 + 0.2e1 * t66) * t44;
t64 = t49 * t30 - t47 * t32;
t46 = sin(qJ(5));
t48 = cos(qJ(5));
t60 = -t44 * mrSges(5,1) + t42 * mrSges(5,2);
t57 = qJDD(3) * mrSges(5,3) + t50 * t60;
t10 = (t44 * t71 - t19 - t67) * t42 + t69;
t11 = -t40 * t71 + t44 * t67 + t65;
t58 = -t42 * t46 + t44 * t48;
t24 = t58 * qJD(3);
t59 = t42 * t48 + t44 * t46;
t25 = t59 * qJD(3);
t17 = -t24 * mrSges(6,1) + t25 * mrSges(6,2);
t21 = t24 * qJD(5) + t59 * qJDD(3);
t22 = -qJD(5) * mrSges(6,2) + t24 * mrSges(6,3);
t8 = m(6) * (t48 * t10 - t46 * t11) - t21 * mrSges(6,3) + qJDD(5) * mrSges(6,1) - t25 * t17 + qJD(5) * t22;
t20 = -t25 * qJD(5) + t58 * qJDD(3);
t23 = qJD(5) * mrSges(6,1) - t25 * mrSges(6,3);
t9 = m(6) * (t46 * t10 + t48 * t11) + t20 * mrSges(6,3) - qJDD(5) * mrSges(6,2) + t24 * t17 - qJD(5) * t23;
t5 = m(5) * t69 + t46 * t9 + t48 * t8 + (-m(5) * t19 - t57) * t42;
t6 = m(5) * t65 + t57 * t44 - t46 * t8 + t48 * t9;
t4 = m(4) * t70 - t50 * mrSges(4,1) - qJDD(3) * mrSges(4,2) - t42 * t5 + t44 * t6;
t61 = qJDD(4) - t64;
t54 = t20 * mrSges(6,1) + t24 * t22 - m(6) * ((-pkin(4) * t44 - pkin(3)) * qJDD(3) + (-t68 * pkin(6) - qJ(4)) * t50 + t61) - t25 * t23 - t21 * mrSges(6,2);
t51 = m(5) * (-qJDD(3) * pkin(3) - t50 * qJ(4) + t61) - t54;
t7 = m(4) * t64 + (-mrSges(4,2) + t72) * t50 + (mrSges(4,1) - t60) * qJDD(3) - t51;
t63 = m(3) * t32 + t49 * t4 - t47 * t7;
t56 = m(3) * t30 + t47 * t4 + t49 * t7;
t55 = m(4) * t41 + t42 * t6 + t44 * t5;
t53 = -m(3) * t41 - t55;
t52 = -m(2) * t41 + t53;
t2 = m(2) * t32 + t63;
t1 = m(2) * t62 - t56;
t3 = [-m(1) * g(1) - t43 * t1 + t45 * t2, t2, t63, t4, t6, t9; -m(1) * g(2) + t45 * t1 + t43 * t2, t1, t53, t7, t5, t8; -m(1) * g(3) + t52, t52, t56, t55, t60 * qJDD(3) - t50 * t72 + t51, -t54;];
f_new = t3;
