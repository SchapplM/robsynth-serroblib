% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PPRPR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
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
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PPRPR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR3_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR3_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR3_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR3_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR3_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:04:57
% EndTime: 2019-12-05 15:04:57
% DurationCPUTime: 0.20s
% Computational Cost: add. (1083->91), mult. (1663->125), div. (0->0), fcn. (1116->10), ass. (0->48)
t211 = sin(pkin(7));
t214 = cos(pkin(7));
t205 = -t214 * g(1) - t211 * g(2);
t208 = -g(3) + qJDD(1);
t210 = sin(pkin(8));
t213 = cos(pkin(8));
t193 = t213 * t205 + t210 * t208;
t204 = -t211 * g(1) + t214 * g(2) + qJDD(2);
t216 = sin(qJ(3));
t218 = cos(qJ(3));
t188 = -t216 * t193 + t218 * t204;
t186 = qJDD(3) * pkin(3) + t188;
t189 = t218 * t193 + t216 * t204;
t219 = qJD(3) ^ 2;
t187 = -t219 * pkin(3) + t189;
t209 = sin(pkin(9));
t212 = cos(pkin(9));
t183 = t209 * t186 + t212 * t187;
t181 = -t219 * pkin(4) + qJDD(3) * pkin(6) + t183;
t222 = t210 * t205 - t213 * t208;
t191 = qJDD(4) + t222;
t215 = sin(qJ(5));
t217 = cos(qJ(5));
t178 = -t215 * t181 + t217 * t191;
t201 = (-t217 * mrSges(6,1) + t215 * mrSges(6,2)) * qJD(3);
t224 = qJD(3) * qJD(5);
t202 = t215 * qJDD(3) + t217 * t224;
t225 = qJD(3) * t217;
t207 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t225;
t226 = qJD(3) * t215;
t176 = m(6) * t178 + qJDD(5) * mrSges(6,1) - t202 * mrSges(6,3) + qJD(5) * t207 - t201 * t226;
t179 = t217 * t181 + t215 * t191;
t203 = t217 * qJDD(3) - t215 * t224;
t206 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t226;
t177 = m(6) * t179 - qJDD(5) * mrSges(6,2) + t203 * mrSges(6,3) - qJD(5) * t206 + t201 * t225;
t223 = -t215 * t176 + t217 * t177;
t172 = m(5) * t183 - t219 * mrSges(5,1) - qJDD(3) * mrSges(5,2) + t223;
t182 = t212 * t186 - t209 * t187;
t180 = -qJDD(3) * pkin(4) - t219 * pkin(6) - t182;
t220 = -m(6) * t180 + t203 * mrSges(6,1) - t202 * mrSges(6,2) - t206 * t226 + t207 * t225;
t174 = m(5) * t182 + qJDD(3) * mrSges(5,1) - t219 * mrSges(5,2) + t220;
t227 = t209 * t172 + t212 * t174;
t221 = m(5) * t191 + t217 * t176 + t215 * t177;
t196 = Ifges(6,5) * qJD(5) + (t215 * Ifges(6,1) + t217 * Ifges(6,4)) * qJD(3);
t195 = Ifges(6,6) * qJD(5) + (t215 * Ifges(6,4) + t217 * Ifges(6,2)) * qJD(3);
t170 = m(4) * t189 - t219 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t212 * t172 - t209 * t174;
t169 = m(4) * t188 + qJDD(3) * mrSges(4,1) - t219 * mrSges(4,2) + t227;
t1 = [m(2) * t208 + t210 * (m(3) * t193 - t216 * t169 + t218 * t170) + t213 * ((-m(3) - m(4)) * t222 - t221); m(3) * t204 + t218 * t169 + t216 * t170; mrSges(4,1) * t188 - mrSges(4,2) * t189 + mrSges(5,1) * t182 - mrSges(5,2) * t183 + t215 * (mrSges(6,2) * t180 - mrSges(6,3) * t178 + Ifges(6,1) * t202 + Ifges(6,4) * t203 + Ifges(6,5) * qJDD(5) - qJD(5) * t195) + t217 * (-mrSges(6,1) * t180 + mrSges(6,3) * t179 + Ifges(6,4) * t202 + Ifges(6,2) * t203 + Ifges(6,6) * qJDD(5) + qJD(5) * t196) + pkin(4) * t220 + pkin(6) * t223 + pkin(3) * t227 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3); t221; mrSges(6,1) * t178 - mrSges(6,2) * t179 + Ifges(6,5) * t202 + Ifges(6,6) * t203 + Ifges(6,3) * qJDD(5) + (t215 * t195 - t217 * t196) * qJD(3);];
tauJ = t1;
