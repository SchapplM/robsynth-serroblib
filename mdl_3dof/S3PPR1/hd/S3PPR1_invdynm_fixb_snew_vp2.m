% Calculate vector of cutting torques with Newton-Euler for
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
% m [3x4]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:21
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S3PPR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PPR1_invdynm_fixb_snew_vp2: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PPR1_invdynm_fixb_snew_vp2: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3PPR1_invdynm_fixb_snew_vp2: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3PPR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PPR1_invdynm_fixb_snew_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3PPR1_invdynm_fixb_snew_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3PPR1_invdynm_fixb_snew_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3PPR1_invdynm_fixb_snew_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:21:17
% EndTime: 2019-05-04 18:21:17
% DurationCPUTime: 0.12s
% Computational Cost: add. (391->58), mult. (438->62), div. (0->0), fcn. (152->2), ass. (0->26)
t42 = -g(2) + qJDD(1);
t43 = -g(1) + qJDD(2);
t44 = sin(qJ(3));
t45 = cos(qJ(3));
t36 = -t44 * t42 + t45 * t43;
t48 = qJD(3) ^ 2;
t32 = m(4) * t36 + qJDD(3) * mrSges(4,1) - t48 * mrSges(4,2);
t37 = t45 * t42 + t44 * t43;
t33 = m(4) * t37 - t48 * mrSges(4,1) - qJDD(3) * mrSges(4,2);
t56 = -t44 * t32 + t45 * t33;
t20 = m(3) * t42 + t56;
t30 = -mrSges(4,1) * g(3) + mrSges(4,3) * t37 + t48 * Ifges(4,5) + Ifges(4,6) * qJDD(3);
t31 = mrSges(4,2) * g(3) - mrSges(4,3) * t36 + Ifges(4,5) * qJDD(3) - t48 * Ifges(4,6);
t59 = pkin(3) * t56 + t45 * t30 + t44 * t31;
t62 = pkin(1) * t20 + (-mrSges(3,2) + mrSges(2,1)) * t42 + t59;
t61 = m(3) + m(4);
t22 = t45 * t32 + t44 * t33;
t58 = mrSges(4,1) * t36 - mrSges(4,2) * t37 + Ifges(4,3) * qJDD(3);
t57 = -pkin(2) * m(4) - mrSges(3,1);
t55 = mrSges(2,3) - t57;
t54 = m(3) * t43 + t22;
t52 = mrSges(3,2) * t43 + mrSges(3,3) * g(3) - pkin(3) * t22 - t44 * t30 + t45 * t31;
t51 = -mrSges(3,1) * t43 + mrSges(3,3) * t42 - pkin(2) * t22 - t58;
t50 = qJ(2) * t61 * g(3) + mrSges(2,1) * g(1) - pkin(1) * t54 + t52;
t49 = -mrSges(2,2) * t42 + qJ(2) * t20 + t51;
t1 = [mrSges(1,3) * g(2) + (qJ(1) * (-m(2) - t61) - mrSges(1,2) - t55) * g(3) + t62, -mrSges(2,3) * g(1) - t49, t52, t31; -mrSges(1,3) * g(1) + (mrSges(1,1) - mrSges(2,2)) * g(3) + t50, g(3) * t55 - t62, t51, t30; -qJ(1) * t54 - mrSges(1,1) * g(2) + (qJ(1) * m(2) + mrSges(1,2) + mrSges(2,3)) * g(1) + t49, -mrSges(2,2) * g(3) + t50, -mrSges(3,2) * t42 + g(3) * t57 + t59, t58;];
m_new  = t1;
