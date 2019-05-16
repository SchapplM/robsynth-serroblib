% Calculate vector of cutting forces with Newton-Euler
% S3PPR1
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% qJDD [3x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d3]';
% m_mdh [4x1]
%   mass of all robot links (including the base)
% mrSges [4x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% f_new [3x4]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:21
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S3PPR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PPR1_invdynf_fixb_snew_vp2: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PPR1_invdynf_fixb_snew_vp2: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3PPR1_invdynf_fixb_snew_vp2: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3PPR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PPR1_invdynf_fixb_snew_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3PPR1_invdynf_fixb_snew_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3PPR1_invdynf_fixb_snew_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3PPR1_invdynf_fixb_snew_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:21:17
% EndTime: 2019-05-04 18:21:17
% DurationCPUTime: 0.05s
% Computational Cost: add. (102->22), mult. (118->25), div. (0->0), fcn. (40->2), ass. (0->13)
t18 = m(2) + m(3);
t10 = -g(1) + qJDD(2);
t11 = sin(qJ(3));
t12 = cos(qJ(3));
t14 = qJD(3) ^ 2;
t9 = -g(2) + qJDD(1);
t4 = m(4) * (t12 * t10 - t11 * t9) + qJDD(3) * mrSges(4,1) - t14 * mrSges(4,2);
t5 = m(4) * (t11 * t10 + t12 * t9) - qJDD(3) * mrSges(4,2) - t14 * mrSges(4,1);
t17 = m(3) * t10 + t11 * t5 + t12 * t4;
t16 = m(3) * t9 - t11 * t4 + t12 * t5;
t15 = m(2) * t9 + t16;
t13 = m(4) * g(3);
t1 = [(-m(1) - m(2)) * g(1) + t17, t18 * g(3) + t13, t16, t5; -m(1) * g(2) + t15, m(2) * g(1) - t17, -m(3) * g(3) - t13, t4; -t13 + (-m(1) - t18) * g(3), t15, t17, t13;];
f_new  = t1;
