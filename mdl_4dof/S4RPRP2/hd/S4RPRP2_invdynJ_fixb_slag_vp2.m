% Calculate vector of inverse dynamics joint torques for
% S4RPRP2
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:50
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau = S4RPRP2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP2_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP2_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP2_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP2_invdynJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP2_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP2_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP2_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:49:23
% EndTime: 2018-11-14 13:49:24
% DurationCPUTime: 0.62s
% Computational Cost: add. (446->103), mult. (706->115), div. (0->0), fcn. (288->4), ass. (0->46)
t40 = -pkin(1) - pkin(2);
t23 = t40 * qJD(1) + qJD(2);
t36 = sin(qJ(3));
t38 = cos(qJ(3));
t48 = qJ(2) * qJD(1);
t14 = t23 * t36 + t38 * t48;
t71 = qJD(1) - qJD(3);
t72 = t14 * t71;
t53 = mrSges(4,2) + mrSges(5,2);
t70 = Ifges(5,3) + Ifges(4,3);
t13 = t38 * t23 - t36 * t48;
t20 = -qJ(2) * t36 + t38 * t40;
t69 = t71 * t36;
t68 = -m(5) - m(3) - m(4);
t22 = t40 * qJDD(1) + qJDD(2);
t47 = qJD(1) * qJD(2);
t24 = qJ(2) * qJDD(1) + t47;
t4 = t13 * qJD(3) + t22 * t36 + t24 * t38;
t37 = sin(qJ(1));
t39 = cos(qJ(1));
t43 = -g(1) * t37 + g(2) * t39;
t67 = t4 * t36 - t38 * t72 + t43;
t15 = -t37 * t36 - t39 * t38;
t16 = t36 * t39 - t37 * t38;
t66 = t16 * g(1) - t15 * g(2);
t64 = t66 - t72;
t62 = m(5) * pkin(3);
t11 = qJD(2) * t38 + t20 * qJD(3);
t21 = qJ(2) * t38 + t36 * t40;
t61 = t14 * t11 + t4 * t21;
t34 = -qJDD(1) + qJDD(3);
t60 = pkin(3) * t34;
t55 = mrSges(2,1) + mrSges(3,1);
t54 = mrSges(4,1) + mrSges(5,1);
t45 = t24 + t47;
t5 = -t14 * qJD(3) + t38 * t22 - t24 * t36;
t2 = t5 + t60;
t44 = -t5 * mrSges(4,1) - t2 * mrSges(5,1);
t42 = -t36 * t62 + mrSges(2,2) - mrSges(3,3);
t41 = qJD(1) ^ 2;
t30 = -pkin(1) * qJDD(1) + qJDD(2);
t28 = pkin(3) * t38 + pkin(2);
t17 = -pkin(3) + t20;
t12 = -qJD(2) * t36 - t21 * qJD(3);
t8 = -pkin(3) * t71 + t13;
t1 = [-t30 * mrSges(3,1) + t53 * t4 + t45 * mrSges(3,3) + m(5) * (t12 * t8 + t17 * t2 + t61) + m(4) * (t12 * t13 + t20 * t5 + t61) + m(3) * (-pkin(1) * t30 + t45 * qJ(2)) - (-t53 * t11 + t54 * t12) * t71 + (mrSges(3,1) * pkin(1) + mrSges(3,3) * qJ(2) + Ifges(3,2) + Ifges(2,3)) * qJDD(1) + ((-m(4) * pkin(2) - m(5) * t28 - t55) * t39 + t42 * t37 + t53 * t16 + t54 * t15 + t68 * (t39 * pkin(1) + t37 * qJ(2))) * g(2) + (-t54 * t16 + t53 * t15 + (-m(5) * (-pkin(1) - t28) + m(3) * pkin(1) - m(4) * t40 + t55) * t37 + (t68 * qJ(2) + t42) * t39) * g(1) + (mrSges(4,1) * t20 + mrSges(5,1) * t17 - t53 * t21 - t70) * t34 + t44; -qJDD(1) * mrSges(3,1) - t41 * mrSges(3,3) + (-t53 * t36 + t54 * t38) * t34 + (-t41 * qJ(2) + t30 + t43) * m(3) + (t2 * t38 + t69 * t8 + t67) * m(5) + (t69 * t13 + t5 * t38 + t67) * m(4) - t71 ^ 2 * (t54 * t36 + t53 * t38); -t44 + t70 * t34 + (t2 + t66) * t62 + t64 * mrSges(4,1) + (t60 + t64) * mrSges(5,1) + t53 * (-t15 * g(1) - t16 * g(2) - t13 * t71 - t4) - m(5) * (t13 - t8) * t14; (g(3) + qJDD(4)) * m(5);];
tau  = t1;
