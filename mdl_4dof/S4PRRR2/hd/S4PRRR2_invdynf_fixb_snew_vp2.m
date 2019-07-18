% Calculate vector of cutting forces with Newton-Euler
% S4PRRR2
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
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
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4PRRR2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(2,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR2_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR2_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_invdynf_fixb_snew_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR2_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR2_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR2_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:27:22
% EndTime: 2019-07-18 13:27:22
% DurationCPUTime: 0.16s
% Computational Cost: add. (1292->46), mult. (1681->58), div. (0->0), fcn. (868->6), ass. (0->36)
t44 = -m(1) - m(2);
t30 = sin(qJ(2));
t33 = cos(qJ(2));
t41 = t33 * g(1) + t30 * g(3);
t13 = qJDD(2) * pkin(1) + t41;
t34 = qJD(2) ^ 2;
t38 = t30 * g(1) - t33 * g(3);
t14 = -t34 * pkin(1) + t38;
t29 = sin(qJ(3));
t32 = cos(qJ(3));
t43 = t29 * t13 + t32 * t14;
t27 = g(2) + qJDD(1);
t18 = m(5) * t27;
t42 = m(4) * t27 + t18;
t26 = qJD(2) + qJD(3);
t25 = qJDD(2) + qJDD(3);
t40 = m(3) * t27 + t42;
t24 = t26 ^ 2;
t28 = sin(qJ(4));
t31 = cos(qJ(4));
t37 = t32 * t13 - t29 * t14;
t17 = qJD(4) + t26;
t15 = t17 ^ 2;
t16 = qJDD(4) + t25;
t8 = t25 * pkin(2) + t37;
t9 = -t24 * pkin(2) + t43;
t6 = m(5) * (-t28 * t9 + t31 * t8) + t16 * mrSges(5,1) - t15 * mrSges(5,2);
t7 = m(5) * (t28 * t8 + t31 * t9) - t16 * mrSges(5,2) - t15 * mrSges(5,1);
t4 = m(4) * t37 + t25 * mrSges(4,1) - t24 * mrSges(4,2) + t28 * t7 + t31 * t6;
t5 = m(4) * t43 - t24 * mrSges(4,1) - t25 * mrSges(4,2) - t28 * t6 + t31 * t7;
t2 = m(3) * t41 + qJDD(2) * mrSges(3,1) - t34 * mrSges(3,2) + t29 * t5 + t32 * t4;
t3 = m(3) * t38 - t34 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t29 * t4 + t32 * t5;
t39 = -t30 * t2 + t33 * t3;
t36 = -t33 * t2 - t30 * t3;
t35 = m(2) * t27 + t40;
t1 = [t44 * g(1) + t36, -m(2) * g(3) + t39, t3, t5, t7; -m(1) * g(2) - t35, m(2) * g(1) - t36, t2, t4, t6; t44 * g(3) + t39, t35, t40, t42, t18;];
f_new  = t1;
