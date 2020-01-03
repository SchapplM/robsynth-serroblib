% Calculate vector of cutting forces with Newton-Euler
% S5RPRPR6
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRPR6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR6_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR6_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR6_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR6_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR6_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR6_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR6_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:17:44
% EndTime: 2019-12-31 18:17:44
% DurationCPUTime: 0.32s
% Computational Cost: add. (3669->83), mult. (4936->106), div. (0->0), fcn. (2392->8), ass. (0->52)
t69 = -pkin(3) - pkin(7);
t39 = qJD(1) + qJD(3);
t43 = sin(qJ(5));
t68 = t39 * t43;
t46 = cos(qJ(5));
t67 = t39 * t46;
t66 = mrSges(4,1) - mrSges(5,2);
t65 = -mrSges(4,2) + mrSges(5,3);
t45 = sin(qJ(1));
t48 = cos(qJ(1));
t61 = t45 * g(1) - t48 * g(2);
t29 = qJDD(1) * pkin(1) + t61;
t49 = qJD(1) ^ 2;
t57 = -t48 * g(1) - t45 * g(2);
t30 = -t49 * pkin(1) + t57;
t41 = sin(pkin(8));
t42 = cos(pkin(8));
t59 = t42 * t29 - t41 * t30;
t18 = qJDD(1) * pkin(2) + t59;
t63 = t41 * t29 + t42 * t30;
t19 = -t49 * pkin(2) + t63;
t44 = sin(qJ(3));
t47 = cos(qJ(3));
t64 = t44 * t18 + t47 * t19;
t62 = qJD(5) * t39;
t60 = t47 * t18 - t44 * t19;
t40 = -g(3) + qJDD(2);
t38 = qJDD(1) + qJDD(3);
t37 = t39 ^ 2;
t51 = -t37 * qJ(4) + qJDD(4) - t60;
t12 = t69 * t38 + t51;
t23 = (mrSges(6,1) * t43 + mrSges(6,2) * t46) * t39;
t25 = t46 * t38 - t43 * t62;
t31 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t68;
t8 = m(6) * (t46 * t12 - t43 * t40) - t25 * mrSges(6,3) + qJDD(5) * mrSges(6,1) - t23 * t67 + qJD(5) * t31;
t24 = -t43 * t38 - t46 * t62;
t32 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t67;
t9 = m(6) * (t43 * t12 + t46 * t40) + t24 * mrSges(6,3) - qJDD(5) * mrSges(6,2) - t23 * t68 - qJD(5) * t32;
t58 = m(5) * t40 - t43 * t8 + t46 * t9;
t56 = m(4) * t40 + t58;
t55 = m(3) * t40 + t56;
t52 = t38 * qJ(4) + 0.2e1 * qJD(4) * t39 + t64;
t54 = -t24 * mrSges(6,1) + m(6) * (t69 * t37 + t52) + t31 * t68 + t32 * t67 + t25 * mrSges(6,2);
t53 = -m(5) * (-t38 * pkin(3) + t51) - t43 * t9 - t46 * t8;
t50 = -m(5) * (t37 * pkin(3) - t52) + t54;
t6 = m(4) * t64 - t66 * t37 + t65 * t38 + t50;
t5 = m(4) * t60 + t65 * t37 + t66 * t38 + t53;
t4 = m(3) * t63 - t49 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t44 * t5 + t47 * t6;
t3 = m(3) * t59 + qJDD(1) * mrSges(3,1) - t49 * mrSges(3,2) + t44 * t6 + t47 * t5;
t2 = m(2) * t57 - t49 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t41 * t3 + t42 * t4;
t1 = m(2) * t61 + qJDD(1) * mrSges(2,1) - t49 * mrSges(2,2) + t42 * t3 + t41 * t4;
t7 = [-m(1) * g(1) - t45 * t1 + t48 * t2, t2, t4, t6, t58, t9; -m(1) * g(2) + t48 * t1 + t45 * t2, t1, t3, t5, -t37 * mrSges(5,2) - t38 * mrSges(5,3) - t50, t8; (-m(1) - m(2)) * g(3) + t55, -m(2) * g(3) + t55, t55, t56, t38 * mrSges(5,2) - t37 * mrSges(5,3) - t53, t54;];
f_new = t7;
