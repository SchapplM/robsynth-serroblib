% Calculate vector of cutting forces with Newton-Euler
% S4PPRR2
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta2]';
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
% Datum: 2019-05-04 18:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4PPRR2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR2_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR2_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR2_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR2_invdynf_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR2_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR2_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR2_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:50:36
% EndTime: 2019-05-04 18:50:36
% DurationCPUTime: 0.14s
% Computational Cost: add. (1054->39), mult. (1309->51), div. (0->0), fcn. (868->6), ass. (0->31)
t39 = -m(1) - m(2);
t25 = -g(2) + qJDD(1);
t26 = sin(pkin(6));
t27 = cos(pkin(6));
t15 = t26 * g(1) + t27 * t25;
t16 = -t27 * g(1) + t26 * t25;
t29 = sin(qJ(3));
t31 = cos(qJ(3));
t38 = t29 * t15 + t31 * t16;
t24 = -g(3) + qJDD(2);
t17 = m(5) * t24;
t37 = m(4) * t24 + t17;
t28 = sin(qJ(4));
t30 = cos(qJ(4));
t32 = qJD(3) ^ 2;
t33 = t31 * t15 - t29 * t16;
t10 = qJDD(3) * pkin(3) + t33;
t11 = -t32 * pkin(3) + t38;
t23 = qJD(3) + qJD(4);
t21 = t23 ^ 2;
t22 = qJDD(3) + qJDD(4);
t8 = m(5) * (t30 * t10 - t28 * t11) + t22 * mrSges(5,1) - t21 * mrSges(5,2);
t9 = m(5) * (t28 * t10 + t30 * t11) - t22 * mrSges(5,2) - t21 * mrSges(5,1);
t6 = m(4) * t33 + qJDD(3) * mrSges(4,1) - t32 * mrSges(4,2) + t28 * t9 + t30 * t8;
t7 = m(4) * t38 - t32 * mrSges(4,1) - qJDD(3) * mrSges(4,2) - t28 * t8 + t30 * t9;
t4 = m(3) * t15 + t29 * t7 + t31 * t6;
t5 = m(3) * t16 - t29 * t6 + t31 * t7;
t36 = m(2) * t25 + t26 * t5 + t27 * t4;
t35 = m(3) * t24 + t37;
t34 = -t26 * t4 + t27 * t5;
t1 = [g(1) * t39 + t34, -m(2) * g(1) + t34, t5, t7, t9; -m(1) * g(2) + t36, m(2) * g(3) - t35, t4, t6, t8; g(3) * t39 + t35, t36, t35, t37, t17;];
f_new  = t1;
