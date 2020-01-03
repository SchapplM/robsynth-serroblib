% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRPP5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP5_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP5_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP5_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP5_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:57:38
% EndTime: 2019-12-31 20:57:42
% DurationCPUTime: 1.75s
% Computational Cost: add. (3430->190), mult. (6898->239), div. (0->0), fcn. (6251->4), ass. (0->110)
t215 = sin(qJ(3));
t216 = sin(qJ(2));
t217 = cos(qJ(3));
t218 = cos(qJ(2));
t121 = t215 * t216 - t217 * t218;
t122 = -t215 * t218 - t216 * t217;
t143 = -pkin(3) - pkin(4);
t138 = -pkin(2) * t218 - pkin(1);
t186 = t122 * qJ(4);
t155 = t138 + t186;
t26 = t143 * t121 - t155;
t241 = m(6) * t26 - mrSges(6,1) * t121 - mrSges(6,2) * t122;
t177 = t218 * pkin(6);
t126 = pkin(7) * t218 + t177;
t174 = t216 * pkin(6);
t154 = -pkin(7) * t216 - t174;
t233 = t217 * t126 + t215 * t154;
t243 = t233 * mrSges(5,1);
t244 = t233 * mrSges(4,1);
t83 = -t215 * t126 + t217 * t154;
t245 = t83 * mrSges(5,3);
t246 = t83 * mrSges(4,2);
t240 = t121 * qJ(5) + t233;
t255 = t240 * mrSges(6,1);
t51 = t122 * qJ(5) - t83;
t262 = t51 * mrSges(6,2);
t265 = t245 - t243 - t244 - t246 - t255 - t262;
t264 = t243 / 0.2e1 + t244 / 0.2e1 - t245 / 0.2e1 + t246 / 0.2e1 + t255 / 0.2e1 + t262 / 0.2e1;
t193 = t121 * mrSges(6,3);
t105 = t193 / 0.2e1;
t195 = t121 * mrSges(5,2);
t107 = -t195 / 0.2e1;
t225 = m(6) / 0.2e1;
t232 = -mrSges(5,2) / 0.2e1 + mrSges(6,3) / 0.2e1;
t261 = t232 * t121 + 0.2e1 * t240 * t225 + t105 + t107;
t239 = Ifges(4,4) - Ifges(5,5) - Ifges(6,4);
t258 = -t51 * mrSges(6,3) + (Ifges(4,1) + Ifges(5,1) + Ifges(6,1) - Ifges(4,2) - Ifges(6,2) - Ifges(5,3)) * t122 + t239 * t121;
t257 = -qJ(4) * t51 + t143 * t240;
t176 = t217 * pkin(2);
t137 = -t176 - pkin(3);
t130 = -pkin(4) + t137;
t173 = t215 * pkin(2);
t133 = t173 + qJ(4);
t256 = t130 * t240 - t133 * t51;
t56 = t121 * pkin(3) + t155;
t242 = m(5) * t56 + mrSges(5,1) * t121 + mrSges(5,3) * t122;
t252 = t121 * t176 + t122 * t173;
t201 = mrSges(6,2) + mrSges(5,3);
t220 = m(5) + m(6);
t129 = qJ(4) * t220 + t201;
t251 = qJD(3) * t129;
t250 = t129 * qJD(4);
t229 = 0.2e1 * m(6);
t238 = -pkin(3) * t233 + qJ(4) * t83;
t231 = m(5) / 0.4e1 + m(6) / 0.4e1;
t230 = -t217 * mrSges(4,2) + (-mrSges(4,1) - mrSges(5,1) - mrSges(6,1)) * t215;
t228 = 0.2e1 * qJD(3);
t227 = m(5) / 0.2e1;
t219 = m(5) * t233;
t214 = t122 * pkin(4);
t7 = (-t241 + t242) * t122;
t197 = qJD(1) * t7;
t67 = t122 * mrSges(6,1) - t121 * mrSges(6,2);
t152 = t138 * (-mrSges(4,1) * t122 - mrSges(4,2) * t121) + t51 * t193 + t26 * t67 + t56 * (-mrSges(5,1) * t122 + mrSges(5,3) * t121);
t158 = t239 * t122;
t175 = t216 * pkin(2);
t110 = qJ(4) * t121;
t65 = -t122 * pkin(3) + t110;
t59 = t175 + t65;
t28 = -t59 + t214;
t1 = -pkin(1) * (mrSges(3,1) * t216 + mrSges(3,2) * t218) + m(4) * t138 * t175 + (-mrSges(4,2) * t175 - t158) * t122 + (mrSges(4,1) * t175 + t258) * t121 + t152 + (-Ifges(3,2) + Ifges(3,1)) * t218 * t216 + (-t216 ^ 2 + t218 ^ 2) * Ifges(3,4) + t242 * t59 + t241 * t28;
t196 = t1 * qJD(1);
t191 = t143 * mrSges(6,3);
t45 = -t65 + t214;
t2 = t258 * t121 - t158 * t122 + t241 * t45 + t242 * t65 + t152;
t190 = t2 * qJD(1);
t12 = -m(6) * (t121 * t240 + t122 * t51) + (-t121 ^ 2 - t122 ^ 2) * mrSges(6,3);
t189 = qJD(1) * t12;
t184 = t133 * t121;
t185 = t130 * t122;
t11 = (t184 / 0.4e1 + t185 / 0.4e1 - t28 / 0.4e1) * t229 - t67;
t188 = t11 * qJD(1);
t183 = t133 * t122;
t182 = t143 * t122;
t15 = (t110 / 0.4e1 + t182 / 0.4e1 - t45 / 0.4e1) * t229 - t67;
t181 = t15 * qJD(1);
t58 = m(6) * t122;
t180 = t58 * qJD(1);
t179 = t201 * t176;
t172 = mrSges(5,2) * t183;
t171 = mrSges(6,3) * t183;
t166 = t217 * t133;
t161 = (Ifges(4,6) - Ifges(5,6)) * t122 + (-Ifges(5,4) - Ifges(4,5)) * t121;
t19 = -t179 + (-m(5) * (t137 * t215 + t166) - m(6) * (t130 * t215 + t166) - t230) * pkin(2);
t144 = (-(t173 - t133) * t83 + (t137 + t176) * t233) * t227 + ((t215 * t51 + t217 * t240) * pkin(2) + t256) * t225 + t130 * t105 + t172 / 0.2e1 - t171 / 0.2e1 + t137 * t107 + t252 * t232 - t264;
t145 = -m(5) * t238 / 0.2e1 - m(6) * t257 / 0.2e1 + pkin(3) * t107 - t121 * t191 / 0.2e1 + t232 * t186 + t264;
t4 = t145 + t144;
t157 = t4 * qJD(1) - t19 * qJD(2);
t85 = 0.4e1 * t231 * t133 + t201;
t156 = qJD(2) * t85;
t153 = qJD(2) * t129 + t251;
t150 = Ifges(6,5) * t121 + Ifges(6,6) * t122 + t161;
t139 = qJ(4) * t176;
t61 = t201 + t220 * t173 / 0.2e1 + t231 * (0.2e1 * t173 + 0.4e1 * qJ(4));
t14 = (t182 + t110 + t45) * t225;
t10 = (t184 + t185 + t28) * t225;
t9 = t219 + t261;
t5 = t233 * t227 + t219 / 0.2e1 + t261;
t3 = t144 - t145 + t150;
t6 = [qJD(2) * t1 + qJD(3) * t2 + qJD(4) * t7 - qJD(5) * t12, t196 + (mrSges(3,2) * t174 - mrSges(3,1) * t177 + Ifges(3,5) * t218 - Ifges(3,6) * t216 + m(4) * (t215 * t83 - t217 * t233) * pkin(2) + t150 + m(6) * t256 + m(5) * (t133 * t83 + t137 * t233) + t130 * t193 + t172 - t171 - t137 * t195 + t252 * mrSges(4,3) + t265) * qJD(2) + t3 * qJD(3) + t5 * qJD(4) + t10 * qJD(5), t190 + t3 * qJD(2) + t9 * qJD(4) + t14 * qJD(5) + (t257 * t225 + t238 * t227) * t228 + (t161 + (Ifges(6,6) + (mrSges(5,2) - mrSges(6,3)) * qJ(4)) * t122 + (pkin(3) * mrSges(5,2) + Ifges(6,5) + t191) * t121 + t265) * qJD(3), qJD(2) * t5 + qJD(3) * t9 + t197, qJD(2) * t10 + qJD(3) * t14 - t189; qJD(3) * t4 + qJD(5) * t11 - t196, -qJD(3) * t19 + qJD(4) * t85, t61 * qJD(4) + ((-pkin(3) * t173 + t139) * t227 + (t143 * t173 + t139) * t225) * t228 + t157 + (t230 * pkin(2) + t179) * qJD(3), qJD(3) * t61 + t156, t188; -qJD(2) * t4 + qJD(5) * t15 - t190, -t157 + t250, t250, t153, t181; qJD(5) * t58 - t197, -t156 - t251, -t153, 0, t180; -qJD(2) * t11 - qJD(3) * t15 - qJD(4) * t58 + t189, -t188, -t181, -t180, 0;];
Cq = t6;
