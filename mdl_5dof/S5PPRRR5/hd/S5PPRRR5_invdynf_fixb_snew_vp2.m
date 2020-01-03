% Calculate vector of cutting forces with Newton-Euler
% S5PPRRR5
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
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1]';
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
% Datum: 2019-12-31 17:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PPRRR5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR5_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR5_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR5_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRR5_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR5_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR5_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR5_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:35:40
% EndTime: 2019-12-31 17:35:40
% DurationCPUTime: 0.29s
% Computational Cost: add. (2969->66), mult. (3821->93), div. (0->0), fcn. (2312->8), ass. (0->47)
t33 = qJDD(3) + qJDD(4);
t38 = sin(qJ(5));
t41 = cos(qJ(5));
t34 = qJD(3) + qJD(4);
t56 = qJD(5) * t34;
t19 = t33 * t38 + t41 * t56;
t20 = t33 * t41 - t38 * t56;
t60 = t34 * t38;
t24 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t60;
t59 = t34 * t41;
t25 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t59;
t32 = t34 ^ 2;
t36 = sin(pkin(8));
t37 = cos(pkin(8));
t52 = g(1) * t36 - g(2) * t37;
t26 = qJDD(2) - t52;
t28 = -g(1) * t37 - g(2) * t36;
t40 = sin(qJ(3));
t43 = cos(qJ(3));
t54 = t43 * t26 - t28 * t40;
t16 = qJDD(3) * pkin(3) + t54;
t44 = qJD(3) ^ 2;
t57 = t40 * t26 + t43 * t28;
t17 = -pkin(3) * t44 + t57;
t39 = sin(qJ(4));
t42 = cos(qJ(4));
t50 = t16 * t42 - t17 * t39;
t61 = (t24 * t38 - t25 * t41) * t34 + m(6) * (-pkin(4) * t33 - pkin(7) * t32 - t50) - t20 * mrSges(6,1) + t19 * mrSges(6,2);
t58 = t39 * t16 + t42 * t17;
t13 = -pkin(4) * t32 + pkin(7) * t33 + t58;
t18 = (-mrSges(6,1) * t41 + mrSges(6,2) * t38) * t34;
t35 = g(3) - qJDD(1);
t10 = m(6) * (-t13 * t38 + t35 * t41) - t19 * mrSges(6,3) + qJDD(5) * mrSges(6,1) - t18 * t60 + qJD(5) * t25;
t11 = m(6) * (t13 * t41 + t35 * t38) + t20 * mrSges(6,3) - qJDD(5) * mrSges(6,2) + t18 * t59 - qJD(5) * t24;
t55 = m(5) * t35 + t41 * t10 + t38 * t11;
t6 = m(5) * t58 - t32 * mrSges(5,1) - t33 * mrSges(5,2) - t38 * t10 + t41 * t11;
t7 = m(5) * t50 + t33 * mrSges(5,1) - t32 * mrSges(5,2) - t61;
t4 = m(4) * t54 + qJDD(3) * mrSges(4,1) - t44 * mrSges(4,2) + t39 * t6 + t42 * t7;
t5 = m(4) * t57 - t44 * mrSges(4,1) - qJDD(3) * mrSges(4,2) - t39 * t7 + t42 * t6;
t53 = m(3) * t28 - t4 * t40 + t43 * t5;
t51 = m(4) * t35 + t55;
t48 = -m(3) * t35 - t51;
t47 = -m(2) * t35 + t48;
t46 = m(3) * t26 + t43 * t4 + t40 * t5;
t2 = m(2) * t28 + t53;
t1 = m(2) * t52 - t46;
t3 = [-m(1) * g(1) - t1 * t36 + t2 * t37, t2, t53, t5, t6, t11; -m(1) * g(2) + t1 * t37 + t2 * t36, t1, t48, t4, t7, t10; -m(1) * g(3) + t47, t47, t46, t51, t55, t61;];
f_new = t3;
