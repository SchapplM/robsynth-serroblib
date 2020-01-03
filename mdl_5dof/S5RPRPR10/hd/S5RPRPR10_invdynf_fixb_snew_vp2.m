% Calculate vector of cutting forces with Newton-Euler
% S5RPRPR10
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRPR10_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR10_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR10_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR10_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR10_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR10_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR10_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:25:57
% EndTime: 2019-12-31 18:25:58
% DurationCPUTime: 0.41s
% Computational Cost: add. (4898->87), mult. (6469->109), div. (0->0), fcn. (2312->8), ass. (0->51)
t34 = -qJDD(1) + qJDD(3);
t41 = sin(qJ(5));
t44 = cos(qJ(5));
t35 = -qJD(1) + qJD(3);
t60 = qJD(5) * t35;
t26 = t41 * t34 + t44 * t60;
t27 = t44 * t34 - t41 * t60;
t65 = t35 * t41;
t28 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t65;
t64 = t35 * t44;
t29 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t64;
t33 = t35 ^ 2;
t47 = qJD(1) ^ 2;
t43 = sin(qJ(1));
t46 = cos(qJ(1));
t55 = -t46 * g(1) - t43 * g(2);
t51 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t55;
t66 = -pkin(1) - pkin(2);
t21 = t66 * t47 + t51;
t57 = t43 * g(1) - t46 * g(2);
t48 = -t47 * qJ(2) + qJDD(2) - t57;
t23 = t66 * qJDD(1) + t48;
t42 = sin(qJ(3));
t45 = cos(qJ(3));
t56 = -t42 * t21 + t45 * t23;
t16 = t34 * pkin(3) + t56;
t61 = t45 * t21 + t42 * t23;
t17 = -t33 * pkin(3) + t61;
t39 = sin(pkin(8));
t40 = cos(pkin(8));
t54 = t40 * t16 - t39 * t17;
t68 = (t41 * t28 - t44 * t29) * t35 + m(6) * (-t34 * pkin(4) - t33 * pkin(7) - t54) - t27 * mrSges(6,1) + t26 * mrSges(6,2);
t67 = -m(3) - m(4);
t63 = mrSges(2,1) + mrSges(3,1);
t62 = t39 * t16 + t40 * t17;
t59 = -m(2) + t67;
t13 = -t33 * pkin(4) + t34 * pkin(7) + t62;
t25 = (-mrSges(6,1) * t44 + mrSges(6,2) * t41) * t35;
t38 = g(3) + qJDD(4);
t10 = m(6) * (-t41 * t13 + t44 * t38) - t26 * mrSges(6,3) + qJDD(5) * mrSges(6,1) - t25 * t65 + qJD(5) * t29;
t11 = m(6) * (t44 * t13 + t41 * t38) + t27 * mrSges(6,3) - qJDD(5) * mrSges(6,2) + t25 * t64 - qJD(5) * t28;
t58 = m(5) * t38 + t44 * t10 + t41 * t11;
t6 = m(5) * t62 - t33 * mrSges(5,1) - t34 * mrSges(5,2) - t41 * t10 + t44 * t11;
t7 = m(5) * t54 + t34 * mrSges(5,1) - t33 * mrSges(5,2) - t68;
t4 = m(4) * t56 + t34 * mrSges(4,1) - t33 * mrSges(4,2) + t39 * t6 + t40 * t7;
t5 = m(4) * t61 - t33 * mrSges(4,1) - t34 * mrSges(4,2) - t39 * t7 + t40 * t6;
t52 = -t42 * t4 + t45 * t5 + m(3) * (-t47 * pkin(1) + t51) + qJDD(1) * mrSges(3,3);
t50 = -m(3) * (-qJDD(1) * pkin(1) + t48) - t45 * t4 - t42 * t5;
t2 = m(2) * t57 + (-mrSges(2,2) + mrSges(3,3)) * t47 + t63 * qJDD(1) + t50;
t1 = m(2) * t55 - qJDD(1) * mrSges(2,2) - t63 * t47 + t52;
t3 = [-m(1) * g(1) + t46 * t1 - t43 * t2, t1, -t47 * mrSges(3,1) + t52, t5, t6, t11; -m(1) * g(2) + t43 * t1 + t46 * t2, t2, t67 * g(3) - t58, t4, t7, t10; (-m(1) + t59) * g(3) - t58, t59 * g(3) - t58, -qJDD(1) * mrSges(3,1) - t47 * mrSges(3,3) - t50, m(4) * g(3) + t58, t58, t68;];
f_new = t3;
