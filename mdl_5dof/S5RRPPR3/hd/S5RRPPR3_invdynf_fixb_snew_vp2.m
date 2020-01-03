% Calculate vector of cutting forces with Newton-Euler
% S5RRPPR3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRPPR3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR3_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR3_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR3_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR3_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR3_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR3_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR3_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:26:30
% EndTime: 2019-12-31 19:26:31
% DurationCPUTime: 0.34s
% Computational Cost: add. (3940->84), mult. (4936->106), div. (0->0), fcn. (2392->8), ass. (0->52)
t68 = -m(2) - m(3);
t67 = -pkin(3) - pkin(7);
t38 = qJD(1) + qJD(2);
t42 = sin(qJ(5));
t66 = t38 * t42;
t45 = cos(qJ(5));
t65 = t38 * t45;
t64 = mrSges(4,1) - mrSges(5,2);
t63 = -mrSges(4,2) + mrSges(5,3);
t37 = qJDD(1) + qJDD(2);
t44 = sin(qJ(1));
t47 = cos(qJ(1));
t59 = t44 * g(1) - t47 * g(2);
t29 = qJDD(1) * pkin(1) + t59;
t48 = qJD(1) ^ 2;
t55 = -t47 * g(1) - t44 * g(2);
t30 = -t48 * pkin(1) + t55;
t43 = sin(qJ(2));
t46 = cos(qJ(2));
t57 = t46 * t29 - t43 * t30;
t18 = t37 * pkin(2) + t57;
t36 = t38 ^ 2;
t61 = t43 * t29 + t46 * t30;
t19 = -t36 * pkin(2) + t61;
t40 = sin(pkin(8));
t41 = cos(pkin(8));
t62 = t40 * t18 + t41 * t19;
t60 = qJD(5) * t38;
t58 = t41 * t18 - t40 * t19;
t39 = -g(3) + qJDD(3);
t50 = -t36 * qJ(4) + qJDD(4) - t58;
t12 = t67 * t37 + t50;
t23 = (mrSges(6,1) * t42 + mrSges(6,2) * t45) * t38;
t25 = t45 * t37 - t42 * t60;
t31 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t66;
t8 = m(6) * (t45 * t12 - t42 * t39) - t25 * mrSges(6,3) + qJDD(5) * mrSges(6,1) - t23 * t65 + qJD(5) * t31;
t24 = -t42 * t37 - t45 * t60;
t32 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t65;
t9 = m(6) * (t42 * t12 + t45 * t39) + t24 * mrSges(6,3) - qJDD(5) * mrSges(6,2) - t23 * t66 - qJD(5) * t32;
t56 = m(5) * t39 - t42 * t8 + t45 * t9;
t54 = m(4) * t39 + t56;
t51 = t37 * qJ(4) + 0.2e1 * qJD(4) * t38 + t62;
t53 = -t24 * mrSges(6,1) + m(6) * (t67 * t36 + t51) + t31 * t66 + t32 * t65 + t25 * mrSges(6,2);
t52 = -m(5) * (-t37 * pkin(3) + t50) - t42 * t9 - t45 * t8;
t49 = -m(5) * (t36 * pkin(3) - t51) + t53;
t6 = m(4) * t62 - t64 * t36 + t63 * t37 + t49;
t5 = m(4) * t58 + t63 * t36 + t64 * t37 + t52;
t4 = m(3) * t61 - t36 * mrSges(3,1) - t37 * mrSges(3,2) - t40 * t5 + t41 * t6;
t3 = m(3) * t57 + t37 * mrSges(3,1) - t36 * mrSges(3,2) + t40 * t6 + t41 * t5;
t2 = m(2) * t55 - t48 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t43 * t3 + t46 * t4;
t1 = m(2) * t59 + qJDD(1) * mrSges(2,1) - t48 * mrSges(2,2) + t46 * t3 + t43 * t4;
t7 = [-m(1) * g(1) - t44 * t1 + t47 * t2, t2, t4, t6, t56, t9; -m(1) * g(2) + t47 * t1 + t44 * t2, t1, t3, t5, -t36 * mrSges(5,2) - t37 * mrSges(5,3) - t49, t8; (-m(1) + t68) * g(3) + t54, t68 * g(3) + t54, -m(3) * g(3) + t54, t54, t37 * mrSges(5,2) - t36 * mrSges(5,3) - t52, t53;];
f_new = t7;
