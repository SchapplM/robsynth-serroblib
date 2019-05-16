% Calculate vector of cutting torques with Newton-Euler for
% S3PRP1
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
%   pkin=[a2,a3,d2]';
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
% Datum: 2019-05-04 18:24
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S3PRP1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRP1_invdynm_fixb_snew_vp2: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PRP1_invdynm_fixb_snew_vp2: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3PRP1_invdynm_fixb_snew_vp2: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3PRP1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PRP1_invdynm_fixb_snew_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3PRP1_invdynm_fixb_snew_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3PRP1_invdynm_fixb_snew_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3PRP1_invdynm_fixb_snew_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:23:33
% EndTime: 2019-05-04 18:23:34
% DurationCPUTime: 0.13s
% Computational Cost: add. (499->71), mult. (651->76), div. (0->0), fcn. (162->2), ass. (0->27)
t75 = m(3) + m(4);
t54 = -g(2) + qJDD(1);
t57 = sin(qJ(2));
t58 = cos(qJ(2));
t45 = -t58 * g(1) + t57 * t54;
t62 = qJD(2) ^ 2;
t40 = -t62 * pkin(2) + qJDD(2) * qJ(3) + (2 * qJD(3) * qJD(2)) + t45;
t74 = mrSges(4,2) * t40;
t73 = -mrSges(3,1) - mrSges(4,1);
t72 = -mrSges(1,2) - mrSges(2,3);
t71 = m(4) * t40 + qJDD(2) * mrSges(4,3);
t32 = m(3) * t45 - qJDD(2) * mrSges(3,2) + t73 * t62 + t71;
t44 = t57 * g(1) + t58 * t54;
t42 = -qJDD(2) * pkin(2) - t62 * qJ(3) + qJDD(3) - t44;
t68 = -m(4) * t42 + qJDD(2) * mrSges(4,1) + t62 * mrSges(4,3);
t33 = m(3) * t44 + qJDD(2) * mrSges(3,1) - t62 * mrSges(3,2) + t68;
t70 = t58 * t32 - t57 * t33;
t69 = mrSges(4,2) * t42 + mrSges(4,3) * g(3) + Ifges(4,4) * qJDD(2) + t62 * Ifges(4,6);
t67 = -mrSges(4,1) * t42 + mrSges(4,3) * t40 + Ifges(4,2) * qJDD(2);
t29 = t74 + mrSges(3,3) * t45 + (Ifges(4,4) + Ifges(3,5)) * t62 + (Ifges(3,6) - Ifges(4,6)) * qJDD(2) + (m(4) * pkin(2) - t73) * g(3);
t30 = -mrSges(3,3) * t44 + Ifges(3,5) * qJDD(2) - t62 * Ifges(3,6) + (m(4) * qJ(3) - mrSges(3,2)) * g(3) + t69;
t66 = mrSges(2,2) * g(1) + pkin(3) * t70 + t58 * t29 + t57 * t30 + (pkin(1) * t75 + mrSges(2,1)) * g(3);
t25 = t57 * t32 + t58 * t33;
t65 = mrSges(2,2) * t54 - pkin(3) * t25 - t57 * t29 + t58 * t30;
t64 = -mrSges(3,2) * t45 + qJ(3) * (-t62 * mrSges(4,1) + t71) + pkin(2) * t68 + mrSges(3,1) * t44 + Ifges(3,3) * qJDD(2) + t67;
t63 = mrSges(2,1) * t54 + pkin(1) * t25 + t64;
t1 = [mrSges(1,3) * g(2) + (qJ(1) * (-m(2) - t75) + t72) * g(3) + t65, -mrSges(2,3) * g(3) + t65, t30, t69; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t66, -mrSges(2,3) * g(1) - t63, t29, t67; -qJ(1) * t70 - mrSges(1,1) * g(2) + (qJ(1) * m(2) - t72) * g(1) + t63, t66, t64, -mrSges(4,1) * g(3) - t62 * Ifges(4,4) + Ifges(4,6) * qJDD(2) - t74;];
m_new  = t1;
