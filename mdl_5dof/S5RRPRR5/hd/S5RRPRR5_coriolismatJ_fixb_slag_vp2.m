% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% m [6x1]
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:03
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRR5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR5_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR5_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR5_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:02:16
% EndTime: 2022-01-20 11:02:25
% DurationCPUTime: 2.89s
% Computational Cost: add. (9988->217), mult. (19245->281), div. (0->0), fcn. (20855->8), ass. (0->135)
t205 = sin(qJ(2));
t293 = t205 * pkin(1);
t194 = qJ(3) + t293;
t201 = sin(pkin(9));
t173 = (-pkin(7) - t194) * t201;
t202 = cos(pkin(9));
t198 = t202 * pkin(7);
t174 = t202 * t194 + t198;
t204 = sin(qJ(4));
t207 = cos(qJ(4));
t136 = t204 * t173 + t207 * t174;
t183 = -t204 * t201 + t207 * t202;
t295 = t183 * pkin(8);
t111 = t136 + t295;
t203 = sin(qJ(5));
t206 = cos(qJ(5));
t135 = t207 * t173 - t204 * t174;
t184 = t207 * t201 + t204 * t202;
t180 = t184 * pkin(8);
t307 = t135 - t180;
t317 = -t203 * t111 + t206 * t307;
t59 = t206 * t111 + t203 * t307;
t338 = -t59 * mrSges(6,1) - t317 * mrSges(6,2);
t155 = t203 * t183 + t206 * t184;
t227 = t206 * t183 - t203 * t184;
t343 = Ifges(6,5) * t227 - Ifges(6,6) * t155;
t12 = t343 + t338;
t345 = t12 * qJD(5);
t188 = (-pkin(7) - qJ(3)) * t201;
t190 = t202 * qJ(3) + t198;
t163 = t204 * t188 + t207 * t190;
t129 = t163 + t295;
t162 = t207 * t188 - t204 * t190;
t306 = t162 - t180;
t318 = -t203 * t129 + t206 * t306;
t93 = t206 * t129 + t203 * t306;
t339 = -t93 * mrSges(6,1) - t318 * mrSges(6,2);
t16 = t343 + t339;
t344 = t16 * qJD(5);
t195 = -t202 * pkin(3) - pkin(2);
t167 = -t183 * pkin(4) + t195;
t313 = t227 * mrSges(6,2);
t321 = t155 * mrSges(6,1);
t326 = t321 + t313;
t331 = t167 * t326;
t208 = cos(qJ(2));
t292 = t208 * pkin(1);
t164 = t167 - t292;
t332 = t164 * t326;
t337 = t331 / 0.2e1 + t332 / 0.2e1;
t294 = t184 * pkin(4);
t335 = m(6) * t294;
t168 = t184 * t292;
t169 = t183 * t292;
t114 = -t206 * t168 - t203 * t169;
t115 = -t203 * t168 + t206 * t169;
t240 = t114 * mrSges(6,1) / 0.2e1 - t115 * mrSges(6,2) / 0.2e1;
t327 = t240 - t168 * mrSges(5,1) / 0.2e1 - t169 * mrSges(5,2) / 0.2e1;
t231 = Ifges(6,4) * t227 ^ 2 + (-Ifges(6,4) * t155 + (Ifges(6,1) - Ifges(6,2)) * t227) * t155;
t323 = -t184 / 0.2e1;
t233 = t313 / 0.2e1;
t320 = t155 * mrSges(6,3);
t238 = t201 ^ 2 + t202 ^ 2;
t309 = t238 * mrSges(4,3);
t319 = -mrSges(3,2) + t309;
t312 = t227 * mrSges(6,3);
t237 = qJD(1) + qJD(2);
t156 = t184 * mrSges(5,1) + t183 * mrSges(5,2);
t253 = t227 * t203;
t257 = t155 * t206;
t301 = m(6) * pkin(4);
t27 = (-t257 / 0.2e1 + t253 / 0.2e1 + t323) * t301 - t326 - t156;
t311 = t237 * t27;
t44 = 0.2e1 * t233 + t321;
t310 = t237 * t44;
t229 = t238 * qJ(3);
t98 = -mrSges(6,1) * t227 + t155 * mrSges(6,2);
t304 = -t202 * mrSges(4,1) - t183 * mrSges(5,1) + t201 * mrSges(4,2) + t184 * mrSges(5,2) - mrSges(3,1) + t98;
t281 = Ifges(5,4) * t184;
t215 = t184 * (Ifges(5,1) * t183 - t281) / 0.2e1 + (Ifges(5,2) * t183 + t281) * t323 + t98 * t294 + (0.2e1 * Ifges(5,4) * t183 + (Ifges(5,1) - Ifges(5,2)) * t184) * t183 / 0.2e1 + t231;
t303 = m(5) / 0.2e1;
t302 = m(6) / 0.2e1;
t300 = t59 / 0.2e1;
t299 = -t317 / 0.2e1;
t296 = pkin(2) * t205;
t283 = -t27 * qJD(4) + t44 * qJD(5);
t45 = t233 - t313 / 0.2e1;
t81 = (t184 + t253 - t257) * t301 / 0.2e1;
t282 = t81 * qJD(4) + t45 * qJD(5);
t279 = pkin(4) * qJD(4);
t187 = t195 - t292;
t130 = t187 * t156;
t41 = t317 * t312;
t42 = t59 * t320;
t3 = t42 - t41 - t215 - t332 + (-t155 * t59 + t227 * t317) * mrSges(6,3) - t164 * t335 - t130;
t265 = t3 * qJD(1);
t7 = -t332 - t231;
t264 = t7 * qJD(1);
t213 = t155 * t320 + t227 * t312 + (t183 ^ 2 + t184 ^ 2) * mrSges(5,3) + t309;
t228 = t238 * t194;
t18 = m(6) * (-t155 * t317 + t227 * t59) + m(5) * (-t135 * t184 + t136 * t183) + m(4) * t228 + t213;
t260 = qJD(1) * t18;
t211 = -t114 * t320 + t115 * t312 + (t168 * t184 + t169 * t183) * mrSges(5,3);
t11 = m(6) * (t114 * t317 + t59 * t115) + m(5) * (-t135 * t168 + t136 * t169) + m(4) * (t228 * t208 - t296) * pkin(1) + t211 + (m(5) * t187 + m(6) * t164 + t304) * t293 + (-m(4) * t293 + t319) * t292;
t259 = t11 * qJD(1);
t245 = t195 * t156;
t234 = t206 * t312;
t230 = (t164 + t167) * t184;
t209 = -t42 / 0.2e1 + t41 / 0.2e1 + (t300 * t155 + t299 * t227) * mrSges(6,3) + t215 + t130 / 0.2e1 + t245 / 0.2e1 + t337;
t217 = m(6) * (t114 * t206 + t115 * t203);
t2 = 0.2e1 * (-t217 / 0.4e1 + m(6) * t230 / 0.4e1) * pkin(4) + t209 - t327;
t6 = t167 * t335 + t215 + t245 + t331;
t226 = t2 * qJD(1) + t6 * qJD(2);
t212 = t231 + t337;
t4 = t212 - t240;
t8 = t231 + t331;
t225 = t4 * qJD(1) + t8 * qJD(2);
t224 = -pkin(4) * t203 * t320 + Ifges(5,5) * t183 - Ifges(5,6) * t184 + t343;
t19 = m(6) * (-t155 * t318 + t227 * t93) + m(5) * (-t162 * t184 + t163 * t183) + m(4) * t229 + t213;
t210 = -m(4) * (t228 + t229) / 0.2e1 - m(5) * ((-t135 - t162) * t184 + (t136 + t163) * t183) / 0.2e1 - m(6) * ((t59 + t93) * t227 - (t317 + t318) * t155) / 0.2e1 - t213;
t219 = (t302 + t303 + m(4) / 0.2e1) * t293;
t9 = t219 + t210;
t223 = qJD(1) * t9 - qJD(2) * t19;
t13 = (t299 + t317 / 0.2e1) * mrSges(6,2) + (-t59 / 0.2e1 + t300) * mrSges(6,1);
t186 = (mrSges(6,1) * t203 + mrSges(6,2) * t206) * pkin(4);
t216 = -t13 * qJD(1) + t186 * qJD(4);
t185 = t186 * qJD(5);
t78 = t81 * qJD(3);
t40 = t44 * qJD(3);
t39 = t45 * qJD(3);
t25 = t27 * qJD(3);
t10 = t219 - t210;
t5 = t212 + t240;
t1 = t209 + t327 + (t217 / 0.2e1 + t230 * t302) * pkin(4);
t14 = [qJD(2) * t11 + qJD(3) * t18 - qJD(4) * t3 - qJD(5) * t7, t10 * qJD(3) + t1 * qJD(4) + t5 * qJD(5) + t259 + (t211 + 0.2e1 * (t114 * t318 + t93 * t115) * t302 + 0.2e1 * (-t162 * t168 + t163 * t169) * t303 + (t319 * t208 + m(4) * (t208 * t229 - t296) + (m(5) * t195 + m(6) * t167 + t304) * t205) * pkin(1)) * qJD(2), qJD(2) * t10 + t260 + t282, -t265 + t1 * qJD(2) + t78 + (-t136 * mrSges(5,1) - t135 * mrSges(5,2) + t224 + t338) * qJD(4) + t345 + (-t234 + m(6) * (t203 * t317 - t206 * t59)) * t279, t5 * qJD(2) + t12 * qJD(4) - t264 + t345 + t39; -qJD(3) * t9 + qJD(4) * t2 + qJD(5) * t4 - t259, qJD(3) * t19 + qJD(4) * t6 + qJD(5) * t8, -t223 + t282, t78 + (-t163 * mrSges(5,1) - t162 * mrSges(5,2) + t224 + t339) * qJD(4) + t344 + (-t234 + m(6) * (t203 * t318 - t206 * t93)) * t279 + t226, t16 * qJD(4) + t225 + t344 + t39; qJD(2) * t9 - t260 + t283, t223 + t283, 0, -t311, t310; -qJD(2) * t2 + qJD(5) * t13 + t25 + t265, -t226 + t25, t311, -t185, -t185 - t216; -qJD(2) * t4 - qJD(4) * t13 + t264 - t40, -t225 - t40, -t310, t216, 0;];
Cq = t14;
