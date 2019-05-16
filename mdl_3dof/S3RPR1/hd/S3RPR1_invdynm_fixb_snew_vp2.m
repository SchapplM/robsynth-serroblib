% Calculate vector of cutting torques with Newton-Euler for
% S3RPR1
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
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
% Datum: 2019-05-04 18:28
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S3RPR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPR1_invdynm_fixb_snew_vp2: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RPR1_invdynm_fixb_snew_vp2: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3RPR1_invdynm_fixb_snew_vp2: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RPR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RPR1_invdynm_fixb_snew_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3RPR1_invdynm_fixb_snew_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3RPR1_invdynm_fixb_snew_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3RPR1_invdynm_fixb_snew_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:28:44
% EndTime: 2019-05-04 18:28:44
% DurationCPUTime: 0.17s
% Computational Cost: add. (1237->90), mult. (1746->102), div. (0->0), fcn. (518->4), ass. (0->39)
t89 = -pkin(1) - pkin(2);
t88 = -m(4) * pkin(2) - mrSges(3,1);
t74 = sin(qJ(1));
t76 = cos(qJ(1));
t57 = t74 * g(1) - t76 * g(2);
t78 = qJD(1) ^ 2;
t58 = -t76 * g(1) - t74 * g(2);
t85 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t58;
t49 = t89 * t78 + t85;
t84 = -t78 * qJ(2) + qJDD(2) - t57;
t52 = t89 * qJDD(1) + t84;
t73 = sin(qJ(3));
t75 = cos(qJ(3));
t47 = -t73 * t49 + t75 * t52;
t64 = -qJD(1) + qJD(3);
t62 = t64 ^ 2;
t63 = -qJDD(1) + qJDD(3);
t44 = m(4) * t47 + t63 * mrSges(4,1) - (t62 * mrSges(4,2));
t48 = t75 * t49 + t73 * t52;
t45 = m(4) * t48 - t62 * mrSges(4,1) - t63 * mrSges(4,2);
t39 = -t73 * t44 + t75 * t45;
t38 = t75 * t44 + t73 * t45;
t53 = -t78 * pkin(1) + t85;
t87 = m(3) * t53 + qJDD(1) * mrSges(3,3) + t39;
t86 = mrSges(4,1) * t47 - mrSges(4,2) * t48 + Ifges(4,3) * t63;
t55 = -qJDD(1) * pkin(1) + t84;
t83 = -m(3) * t55 + qJDD(1) * mrSges(3,1) + t78 * mrSges(3,3) - t38;
t42 = -mrSges(4,1) * g(3) + mrSges(4,3) * t48 + (t62 * Ifges(4,5)) + Ifges(4,6) * t63;
t43 = mrSges(4,2) * g(3) - mrSges(4,3) * t47 + Ifges(4,5) * t63 - t62 * Ifges(4,6);
t82 = mrSges(3,2) * t55 + mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) + t78 * Ifges(3,6) - pkin(4) * t38 - t73 * t42 + t75 * t43;
t81 = -mrSges(3,2) * t53 + pkin(4) * t39 + t75 * t42 + t73 * t43;
t80 = -mrSges(3,1) * t55 + mrSges(3,3) * t53 + Ifges(3,2) * qJDD(1) - pkin(2) * t38 - t86;
t79 = -mrSges(2,2) * t58 + qJ(2) * (-t78 * mrSges(3,1) + t87) + pkin(1) * t83 + mrSges(2,1) * t57 + Ifges(2,3) * qJDD(1) + t80;
t59 = (-m(3) - m(4)) * g(3);
t35 = m(2) * t57 + qJDD(1) * mrSges(2,1) - t78 * mrSges(2,2) + t83;
t34 = m(2) * t58 - qJDD(1) * mrSges(2,2) + (-mrSges(2,1) - mrSges(3,1)) * t78 + t87;
t33 = -mrSges(2,2) * g(3) - mrSges(2,3) * t57 + Ifges(2,5) * qJDD(1) - t78 * Ifges(2,6) - qJ(2) * t59 + t82;
t32 = mrSges(2,3) * t58 - pkin(1) * t59 + (Ifges(3,4) + Ifges(2,5)) * t78 + (Ifges(2,6) - Ifges(3,6)) * qJDD(1) + (mrSges(2,1) - t88) * g(3) - t81;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t76 * t33 - t74 * t32 - pkin(3) * (t74 * t34 + t76 * t35), t33, t82, t43; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t74 * t33 + t76 * t32 + pkin(3) * (t76 * t34 - t74 * t35), t32, t80, t42; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t79, t79, -t78 * Ifges(3,4) + Ifges(3,6) * qJDD(1) + t88 * g(3) + t81, t86;];
m_new  = t1;
