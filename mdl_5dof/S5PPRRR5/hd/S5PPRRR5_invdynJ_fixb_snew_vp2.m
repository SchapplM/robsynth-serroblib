% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PPRRR5
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1]';
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
% Datum: 2019-12-31 17:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PPRRR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR5_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR5_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR5_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRR5_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR5_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR5_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR5_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:35:42
% EndTime: 2019-12-31 17:35:42
% DurationCPUTime: 0.19s
% Computational Cost: add. (984->91), mult. (1255->122), div. (0->0), fcn. (762->8), ass. (0->45)
t209 = qJD(3) + qJD(4);
t213 = sin(qJ(5));
t226 = t209 * t213;
t216 = cos(qJ(5));
t225 = t209 * t216;
t211 = sin(pkin(8));
t212 = cos(pkin(8));
t204 = -t211 * g(1) + t212 * g(2) + qJDD(2);
t205 = -t212 * g(1) - t211 * g(2);
t215 = sin(qJ(3));
t218 = cos(qJ(3));
t189 = t218 * t204 - t215 * t205;
t187 = qJDD(3) * pkin(3) + t189;
t190 = t215 * t204 + t218 * t205;
t219 = qJD(3) ^ 2;
t188 = -t219 * pkin(3) + t190;
t214 = sin(qJ(4));
t217 = cos(qJ(4));
t184 = t214 * t187 + t217 * t188;
t207 = t209 ^ 2;
t208 = qJDD(3) + qJDD(4);
t181 = -t207 * pkin(4) + t208 * pkin(7) + t184;
t210 = g(3) - qJDD(1);
t178 = -t213 * t181 + t216 * t210;
t196 = (-mrSges(6,1) * t216 + mrSges(6,2) * t213) * t209;
t223 = qJD(5) * t209;
t197 = t213 * t208 + t216 * t223;
t203 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t225;
t176 = m(6) * t178 + qJDD(5) * mrSges(6,1) - t197 * mrSges(6,3) + qJD(5) * t203 - t196 * t226;
t179 = t216 * t181 + t213 * t210;
t198 = t216 * t208 - t213 * t223;
t202 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t226;
t177 = m(6) * t179 - qJDD(5) * mrSges(6,2) + t198 * mrSges(6,3) - qJD(5) * t202 + t196 * t225;
t222 = -t213 * t176 + t216 * t177;
t168 = m(5) * t184 - t207 * mrSges(5,1) - t208 * mrSges(5,2) + t222;
t183 = t217 * t187 - t214 * t188;
t180 = -t208 * pkin(4) - t207 * pkin(7) - t183;
t220 = -m(6) * t180 + t198 * mrSges(6,1) - t197 * mrSges(6,2) - t202 * t226 + t203 * t225;
t173 = m(5) * t183 + t208 * mrSges(5,1) - t207 * mrSges(5,2) + t220;
t224 = t214 * t168 + t217 * t173;
t191 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t213 + Ifges(6,6) * t216) * t209;
t192 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t213 + Ifges(6,2) * t216) * t209;
t193 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t213 + Ifges(6,4) * t216) * t209;
t221 = -mrSges(5,2) * t184 + pkin(7) * t222 + t213 * (mrSges(6,2) * t180 - mrSges(6,3) * t178 + Ifges(6,1) * t197 + Ifges(6,4) * t198 + Ifges(6,5) * qJDD(5) - qJD(5) * t192 + t191 * t225) + t216 * (-mrSges(6,1) * t180 + mrSges(6,3) * t179 + Ifges(6,4) * t197 + Ifges(6,2) * t198 + Ifges(6,6) * qJDD(5) + qJD(5) * t193 - t191 * t226) + pkin(4) * t220 + mrSges(5,1) * t183 + Ifges(5,3) * t208;
t1 = [-t216 * t176 - t213 * t177 + (-m(2) - m(3) - m(4) - m(5)) * t210; m(3) * t204 + t215 * (m(4) * t190 - t219 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t217 * t168 - t214 * t173) + t218 * (m(4) * t189 + qJDD(3) * mrSges(4,1) - t219 * mrSges(4,2) + t224); mrSges(4,1) * t189 - mrSges(4,2) * t190 + Ifges(4,3) * qJDD(3) + pkin(3) * t224 + t221; t221; mrSges(6,1) * t178 - mrSges(6,2) * t179 + Ifges(6,5) * t197 + Ifges(6,6) * t198 + Ifges(6,3) * qJDD(5) + (t192 * t213 - t193 * t216) * t209;];
tauJ = t1;
