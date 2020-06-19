% Calculate vector of cutting forces with Newton-Euler
% S2RR3
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% qJDD [2x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,d1,d2]';
% m [3x1]
%   mass of all robot links (including the base)
% mrSges [3x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [3x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% f_new [3x3]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 09:14
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S2RR3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(2,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR3_invdynf_fixb_snew_vp2: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR3_invdynf_fixb_snew_vp2: qJD has to be [2x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [2 1]), ...
  'S2RR3_invdynf_fixb_snew_vp2: qJDD has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S2RR3_invdynf_fixb_snew_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'S2RR3_invdynf_fixb_snew_vp2: m has to be [3x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [3,3]), ...
  'S2RR3_invdynf_fixb_snew_vp2: mrSges has to be [3x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [3 6]), ...
  'S2RR3_invdynf_fixb_snew_vp2: Ifges has to be [3x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-19 09:14:24
% EndTime: 2020-06-19 09:14:25
% DurationCPUTime: 0.14s
% Computational Cost: add. (164->27), mult. (242->37), div. (0->0), fcn. (112->4), ass. (0->18)
t18 = -m(2) - m(3);
t12 = sin(qJ(1));
t14 = cos(qJ(1));
t17 = t12 * g(1) - t14 * g(2);
t16 = -t14 * g(1) - t12 * g(2);
t15 = qJD(1) ^ 2;
t13 = cos(qJ(2));
t11 = sin(qJ(2));
t10 = qJD(1) + qJD(2);
t9 = qJDD(1) + qJDD(2);
t8 = t10 ^ 2;
t6 = -t15 * pkin(1) + t16;
t5 = qJDD(1) * pkin(1) + t17;
t4 = m(3) * (t11 * t5 + t13 * t6) - t9 * mrSges(3,2) - t8 * mrSges(3,1);
t3 = m(3) * (-t11 * t6 + t13 * t5) + t9 * mrSges(3,1) - t8 * mrSges(3,2);
t2 = m(2) * t16 - t15 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t11 * t3 + t13 * t4;
t1 = m(2) * t17 + qJDD(1) * mrSges(2,1) - t15 * mrSges(2,2) + t11 * t4 + t13 * t3;
t7 = [-m(1) * g(1) - t12 * t1 + t14 * t2, t2, t4; -m(1) * g(2) + t14 * t1 + t12 * t2, t1, t3; (-m(1) + t18) * g(3), t18 * g(3), -m(3) * g(3);];
f_new = t7;
