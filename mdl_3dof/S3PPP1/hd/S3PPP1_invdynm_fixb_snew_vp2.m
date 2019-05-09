% Calculate vector of cutting torques with Newton-Euler for
% S3PPP1
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
%   pkin=[a2,a3,theta1]';
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
% m [3x4]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:20
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S3PPP1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PPP1_invdynm_fixb_snew_vp2: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PPP1_invdynm_fixb_snew_vp2: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3PPP1_invdynm_fixb_snew_vp2: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3PPP1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PPP1_invdynm_fixb_snew_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3PPP1_invdynm_fixb_snew_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3PPP1_invdynm_fixb_snew_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3PPP1_invdynm_fixb_snew_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:20:08
% EndTime: 2019-05-04 18:20:08
% DurationCPUTime: 0.09s
% Computational Cost: add. (268->53), mult. (300->54), div. (0->0), fcn. (154->2), ass. (0->26)
t48 = -m(4) * qJ(3) + mrSges(3,2);
t36 = sin(pkin(3));
t37 = cos(pkin(3));
t32 = t36 * g(1) - t37 * g(2);
t29 = qJDD(2) - t32;
t26 = mrSges(4,1) * t29;
t46 = m(4) * pkin(2);
t47 = (mrSges(3,1) + t46) * t29 + t26;
t45 = m(3) + m(4);
t43 = mrSges(4,2) + mrSges(3,3);
t42 = t45 * t29;
t33 = t37 * g(1) + t36 * g(2);
t30 = qJDD(3) - t33;
t35 = -g(3) + qJDD(1);
t34 = mrSges(4,3) * t35;
t40 = -mrSges(4,1) * t30 + t34;
t24 = mrSges(4,2) * t30;
t39 = t24 + (-mrSges(4,3) + t48) * t29;
t27 = m(4) * t30;
t38 = -pkin(1) * t42 + qJ(2) * (-m(3) * t33 + t27) + mrSges(2,1) * t32 + (mrSges(2,2) - mrSges(3,3)) * t33 + t39;
t31 = t45 * t35;
t21 = t27 + (-m(2) - m(3)) * t33;
t20 = m(2) * t32 - t42;
t19 = -mrSges(2,3) * t32 - qJ(2) * t31 + (mrSges(2,2) - t43) * t35 + t47;
t18 = -pkin(1) * t31 - t34 + (-mrSges(3,1) - mrSges(2,3)) * t33 + (mrSges(4,1) + t46) * t30 + (-mrSges(2,1) + t48) * t35;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t37 * t19 - t36 * t18 - qJ(1) * (t37 * t20 + t36 * t21), t19, -mrSges(3,3) * t33 + t39, -mrSges(4,3) * t29 + t24; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t36 * t19 + t37 * t18 + qJ(1) * (-t36 * t20 + t37 * t21), t18, t43 * t35 - t47, t40; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t38, t38, mrSges(3,1) * t33 - mrSges(3,2) * t35 + (-pkin(2) * t30 + qJ(3) * t35) * m(4) + t40, -mrSges(4,2) * t35 + t26;];
m_new  = t1;
