% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRRP5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP5_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP5_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP5_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP5_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:47:49
% EndTime: 2019-12-05 16:47:54
% DurationCPUTime: 1.89s
% Computational Cost: add. (3099->189), mult. (7248->267), div. (0->0), fcn. (7052->6), ass. (0->113)
t131 = sin(qJ(4));
t132 = sin(qJ(3));
t134 = cos(qJ(3));
t203 = cos(qJ(4));
t152 = t203 * t134;
t113 = -t131 * t132 + t152;
t114 = -t131 * t134 - t203 * t132;
t69 = -mrSges(6,1) * t113 - mrSges(6,2) * t114;
t126 = -pkin(3) * t134 - pkin(2);
t90 = -t113 * pkin(4) + t126;
t249 = m(6) * t90 + t69;
t148 = (Ifges(5,6) + Ifges(6,6)) * t114 + (Ifges(5,5) + Ifges(6,5)) * t113;
t212 = -pkin(7) - pkin(6);
t118 = t212 * t132;
t119 = t212 * t134;
t223 = t203 * t118 + t131 * t119;
t231 = t223 * mrSges(5,2);
t230 = t114 * qJ(5) + t223;
t235 = t230 * mrSges(6,2);
t78 = t131 * t118 - t203 * t119;
t238 = t78 * mrSges(5,1);
t58 = t113 * qJ(5) + t78;
t244 = t58 * mrSges(6,1);
t248 = t148 - t231 - t235 - t238 - t244;
t227 = -t78 * mrSges(5,3) - t58 * mrSges(6,3);
t236 = Ifges(5,4) + Ifges(6,4);
t247 = -t236 * t114 + t227;
t246 = -t231 / 0.2e1 - t235 / 0.2e1 - t238 / 0.2e1 - t244 / 0.2e1;
t201 = pkin(3) * t132;
t239 = m(5) * t201;
t232 = mrSges(5,2) + mrSges(6,2);
t160 = t203 * pkin(3);
t125 = t160 + pkin(4);
t147 = t160 - t125;
t233 = mrSges(5,1) + mrSges(6,1);
t229 = t236 * t113;
t226 = Ifges(5,1) + Ifges(6,1);
t225 = Ifges(5,2) + Ifges(6,2);
t222 = t232 * t203;
t178 = t134 * mrSges(4,1);
t220 = t132 * mrSges(4,2) - t178;
t218 = t225 - t226;
t185 = mrSges(6,3) * t113;
t214 = m(6) / 0.2e1;
t217 = (-t185 / 0.2e1 - t58 * t214) * pkin(4) + t246;
t133 = sin(qJ(2));
t100 = t114 * t133;
t169 = t132 * t133;
t99 = t131 * t169 - t133 * t152;
t142 = -t232 * t100 + t233 * t99;
t216 = m(5) / 0.2e1;
t215 = -m(6) / 0.2e1;
t213 = m(6) * pkin(4);
t210 = m(6) * t99;
t135 = cos(qJ(2));
t204 = -t135 / 0.2e1;
t202 = pkin(3) * t131;
t200 = pkin(4) * t114;
t182 = t113 * t99;
t180 = t132 * mrSges(4,1);
t177 = t134 * mrSges(4,2);
t176 = t100 * t114;
t172 = t125 * t114;
t101 = t114 * t135;
t102 = t113 * t135;
t165 = t132 ^ 2 + t134 ^ 2;
t168 = t133 * t135;
t13 = 0.4e1 * (m(5) / 0.4e1 + m(6) / 0.4e1) * (t100 * t101 - t99 * t102 - t168) + m(4) * (-0.1e1 + t165) * t168;
t171 = t13 * qJD(1);
t162 = t113 * t202;
t105 = t113 * mrSges(6,2);
t67 = -t114 * mrSges(6,1) + t105;
t98 = -t200 + t201;
t23 = 0.2e1 * (t162 / 0.4e1 + t172 / 0.4e1 - t98 / 0.4e1) * m(6) - t67;
t167 = t23 * qJD(2);
t158 = mrSges(6,1) + t213;
t38 = t158 * t114 - t105;
t166 = t38 * qJD(2);
t164 = t213 / 0.2e1;
t163 = t114 * t202;
t159 = t125 * t185;
t68 = -mrSges(5,1) * t114 + mrSges(5,2) * t113;
t151 = t227 * t114 - t126 * t68 - t90 * t67;
t150 = t233 * t101 / 0.2e1 - t232 * t102 / 0.2e1;
t149 = t113 * t160;
t70 = -mrSges(5,1) * t113 - mrSges(5,2) * t114;
t1 = (-pkin(2) * mrSges(4,1) - Ifges(4,4) * t132 + pkin(3) * t70) * t132 + t126 * t239 + (-mrSges(4,2) * pkin(2) + Ifges(4,4) * t134 + (Ifges(4,1) - Ifges(4,2)) * t132) * t134 - t151 + t229 * t113 + ((-t226 / 0.2e1 + t225 / 0.2e1 + t218 / 0.2e1) * t113 + t247) * t114 + t249 * t98;
t137 = (t239 / 0.2e1 - t98 * t215) * t135;
t82 = t102 * t202;
t138 = (t101 * t160 + t82) * t216 + (t125 * t101 + t82) * t214 + t150;
t3 = (t68 / 0.2e1 + t67 / 0.2e1) * t135 + t137 + t138;
t146 = -t3 * qJD(1) + t1 * qJD(2);
t2 = -t230 * t185 + (mrSges(6,3) * t230 - t229) * t113 + (t249 * pkin(4) - t218 * t113 - t247) * t114 + t151;
t144 = (t67 + t68) * t204;
t139 = t135 * t200 * t214 + t144;
t6 = (mrSges(5,2) / 0.2e1 + mrSges(6,2) / 0.2e1) * t102 + (-mrSges(5,1) / 0.2e1 - mrSges(6,1) / 0.2e1 - t213 / 0.2e1) * t101 + t139;
t145 = t6 * qJD(1) - t2 * qJD(2);
t20 = m(6) * (t113 * t58 + t114 * t230) + (t113 ^ 2 + t114 ^ 2) * mrSges(6,3);
t27 = (t176 / 0.2e1 - t182 / 0.2e1 - t133 / 0.2e1) * m(6);
t143 = qJD(1) * t27 + qJD(2) * t20;
t21 = (pkin(4) / 0.2e1 + t160 / 0.2e1 - t125 / 0.2e1) * t210;
t48 = t222 * pkin(3) + (-m(6) * t147 + t233) * t202;
t136 = t159 / 0.2e1 - mrSges(6,3) * t149 / 0.2e1 + t147 * t215 * t58 - t246;
t8 = t136 + t217;
t141 = t21 * qJD(1) + t8 * qJD(2) + t48 * qJD(3);
t81 = t100 * t202;
t30 = (t162 + t172 + t98) * t214;
t26 = (t176 - t182 + t133) * t214;
t12 = -t147 * t210 / 0.2e1 + t99 * t164 + t142;
t7 = t101 * t164 + t139 + t150;
t5 = -t136 + t148 + t217;
t4 = (t177 + t180) * t204 + (-t177 / 0.2e1 - t180 / 0.2e1) * t135 - t137 + t138 + t144;
t9 = [qJD(2) * t13, t4 * qJD(3) + t7 * qJD(4) + t26 * qJD(5) + t171 + ((t165 * mrSges(4,3) - mrSges(3,2)) * t135 + (-mrSges(3,1) + t69 + t70 + t220) * t133 + m(4) * (t165 * t135 * pkin(6) - t133 * pkin(2)) + 0.2e1 * (t101 * t223 + t102 * t78 + t126 * t133) * t216 + 0.2e1 * (t101 * t230 + t102 * t58 + t133 * t90) * t214 + (mrSges(5,3) + mrSges(6,3)) * (t101 * t114 + t102 * t113)) * qJD(2), t4 * qJD(2) + (mrSges(4,2) * t169 - t133 * t178 + t142) * qJD(3) + t12 * qJD(4) + 0.2e1 * ((t99 * t160 + t81) * t216 + (t125 * t99 + t81) * t214) * qJD(3), t7 * qJD(2) + t12 * qJD(3) + (pkin(4) * t210 + t142) * qJD(4), t26 * qJD(2); -qJD(3) * t3 + qJD(4) * t6 + qJD(5) * t27 - t171, qJD(3) * t1 - qJD(4) * t2 + qJD(5) * t20, (-t159 + mrSges(6,3) * t163 + m(6) * (-t125 * t58 + t202 * t230) + m(5) * (t131 * t223 - t203 * t78) * pkin(3) + Ifges(4,5) * t134 - Ifges(4,6) * t132 + t220 * pkin(6) + (t163 - t149) * mrSges(5,3) + t248) * qJD(3) + t5 * qJD(4) + t30 * qJD(5) + t146, t5 * qJD(3) + ((-m(6) * t58 - t185) * pkin(4) + t248) * qJD(4) + t145, qJD(3) * t30 + t143; qJD(2) * t3 - qJD(4) * t21, -qJD(4) * t8 + qJD(5) * t23 - t146, -t48 * qJD(4), ((-mrSges(5,1) - t158) * t131 - t222) * qJD(4) * pkin(3) - t141, t167; -qJD(2) * t6 + qJD(3) * t21, qJD(3) * t8 + qJD(5) * t38 - t145, t141, 0, t166; -t27 * qJD(2), -qJD(3) * t23 - qJD(4) * t38 - t143, -t167, -t166, 0;];
Cq = t9;
