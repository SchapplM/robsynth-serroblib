% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PPRPR5
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1]';
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
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PPRPR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR5_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR5_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR5_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR5_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR5_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:33:27
% EndTime: 2019-12-31 17:33:27
% DurationCPUTime: 0.14s
% Computational Cost: add. (407->85), mult. (642->108), div. (0->0), fcn. (322->6), ass. (0->38)
t211 = -pkin(3) - pkin(6);
t196 = sin(pkin(7));
t197 = cos(pkin(7));
t188 = -t196 * g(1) + t197 * g(2) + qJDD(2);
t189 = -t197 * g(1) - t196 * g(2);
t199 = sin(qJ(3));
t201 = cos(qJ(3));
t175 = t199 * t188 + t201 * t189;
t198 = sin(qJ(5));
t210 = qJD(3) * t198;
t200 = cos(qJ(5));
t209 = qJD(3) * t200;
t208 = qJD(3) * qJD(5);
t174 = t201 * t188 - t199 * t189;
t202 = qJD(3) ^ 2;
t206 = -t202 * qJ(4) + qJDD(4) - t174;
t171 = t211 * qJDD(3) + t206;
t193 = g(3) - qJDD(1);
t167 = t200 * t171 - t198 * t193;
t185 = (t198 * mrSges(6,1) + t200 * mrSges(6,2)) * qJD(3);
t187 = t200 * qJDD(3) - t198 * t208;
t190 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t210;
t165 = m(6) * t167 + qJDD(5) * mrSges(6,1) - t187 * mrSges(6,3) + qJD(5) * t190 - t185 * t209;
t168 = t198 * t171 + t200 * t193;
t186 = -t198 * qJDD(3) - t200 * t208;
t191 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t209;
t166 = m(6) * t168 - qJDD(5) * mrSges(6,2) + t186 * mrSges(6,3) - qJD(5) * t191 - t185 * t210;
t207 = t200 * t165 + t198 * t166;
t205 = qJDD(3) * qJ(4) + 0.2e1 * qJD(4) * qJD(3) + t175;
t173 = -qJDD(3) * pkin(3) + t206;
t204 = -m(5) * t173 + t202 * mrSges(5,3) - t207;
t170 = t211 * t202 + t205;
t172 = t202 * pkin(3) - t205;
t203 = -m(5) * t172 + m(6) * t170 - t186 * mrSges(6,1) + t202 * mrSges(5,2) + t187 * mrSges(6,2) + qJDD(3) * mrSges(5,3) + t190 * t210 + t191 * t209;
t178 = Ifges(6,5) * qJD(5) + (t200 * Ifges(6,1) - t198 * Ifges(6,4)) * qJD(3);
t177 = Ifges(6,6) * qJD(5) + (t200 * Ifges(6,4) - t198 * Ifges(6,2)) * qJD(3);
t164 = qJDD(3) * mrSges(5,2) - t204;
t1 = [t198 * t165 - t200 * t166 + (-m(2) - m(3) - m(4) - m(5)) * t193; m(3) * t188 + t199 * (m(4) * t175 - t202 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t203) + t201 * (m(4) * t174 - t202 * mrSges(4,2) + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + t204); mrSges(4,1) * t174 - mrSges(4,2) * t175 + mrSges(5,2) * t173 - mrSges(5,3) * t172 + t200 * (mrSges(6,2) * t170 - mrSges(6,3) * t167 + Ifges(6,1) * t187 + Ifges(6,4) * t186 + Ifges(6,5) * qJDD(5) - qJD(5) * t177) - t198 * (-mrSges(6,1) * t170 + mrSges(6,3) * t168 + Ifges(6,4) * t187 + Ifges(6,2) * t186 + Ifges(6,6) * qJDD(5) + qJD(5) * t178) - pkin(6) * t207 - pkin(3) * t164 + qJ(4) * t203 + (Ifges(4,3) + Ifges(5,1)) * qJDD(3); t164; mrSges(6,1) * t167 - mrSges(6,2) * t168 + Ifges(6,5) * t187 + Ifges(6,6) * t186 + Ifges(6,3) * qJDD(5) + (t200 * t177 + t198 * t178) * qJD(3);];
tauJ = t1;
