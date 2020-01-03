% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PPRRP4_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP4_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP4_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRP4_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRP4_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:34:28
% EndTime: 2019-12-31 17:34:33
% DurationCPUTime: 3.07s
% Computational Cost: add. (4750->192), mult. (10206->286), div. (0->0), fcn. (12192->6), ass. (0->137)
t277 = Icges(5,5) + Icges(6,5);
t280 = -Icges(5,3) - Icges(6,3);
t155 = sin(qJ(4));
t156 = cos(qJ(4));
t169 = Icges(6,5) * t156 - Icges(6,6) * t155;
t171 = Icges(5,5) * t156 - Icges(5,6) * t155;
t279 = t169 + t171;
t204 = Icges(6,4) * t155;
t177 = Icges(6,1) * t156 - t204;
t206 = Icges(5,4) * t155;
t179 = Icges(5,1) * t156 - t206;
t278 = t177 + t179;
t276 = Icges(5,6) + Icges(6,6);
t203 = Icges(6,4) * t156;
t173 = -Icges(6,2) * t155 + t203;
t205 = Icges(5,4) * t156;
t175 = -Icges(5,2) * t155 + t205;
t275 = t173 + t175;
t207 = sin(pkin(7));
t208 = cos(pkin(7));
t230 = sin(qJ(3));
t231 = cos(qJ(3));
t131 = -t207 * t230 - t208 * t231;
t132 = -t207 * t231 + t208 * t230;
t261 = t278 * t131 - t277 * t132;
t262 = t275 * t131 - t276 * t132;
t271 = t279 * t131 + t280 * t132;
t274 = -t271 * t132 + (-t262 * t155 + t261 * t156) * t131;
t260 = t277 * t131 + t278 * t132;
t268 = rSges(6,3) + qJ(5) + pkin(6);
t251 = t268 * t132;
t273 = (t280 * t131 - t279 * t132) * t132;
t270 = -t276 * t131 - t275 * t132;
t269 = t270 * t155 + t271;
t195 = t131 * t156;
t267 = t195 * t260 + t273;
t193 = t132 * t156;
t194 = t132 * t155;
t266 = t271 * t131 + t261 * t193 - t262 * t194;
t196 = t131 * t155;
t265 = t270 * t196 + t267;
t212 = t155 * rSges(6,2);
t215 = rSges(6,1) * t156;
t148 = t212 - t215;
t228 = pkin(4) * t156;
t153 = pkin(3) + t228;
t258 = -t148 + t153;
t257 = t277 * t155 + t276 * t156;
t255 = -t270 * t194 + t274;
t234 = t131 / 0.2e1;
t233 = -t132 / 0.2e1;
t214 = rSges(6,2) * t156;
t164 = (rSges(6,1) + pkin(4)) * t155 + t214;
t161 = t164 * t132;
t72 = t161 * t132;
t252 = t161 * t207;
t172 = Icges(6,2) * t156 + t204;
t174 = Icges(5,2) * t156 + t206;
t250 = -t172 - t174;
t176 = Icges(6,1) * t155 + t203;
t178 = Icges(5,1) * t155 + t205;
t249 = -t176 - t178;
t190 = qJD(4) * t132;
t182 = t155 * rSges(5,1) + rSges(5,2) * t156;
t114 = t182 * t132;
t115 = t182 * t131;
t180 = rSges(6,2) * t194 - t131 * t268;
t185 = -t153 - t215;
t50 = t132 * t185 + t180;
t51 = t258 * t131 - t251;
t128 = t131 * pkin(6);
t184 = rSges(5,2) * t194 - t131 * rSges(5,3);
t216 = rSges(5,1) * t156;
t54 = -t128 + (-pkin(3) - t216) * t132 + t184;
t127 = t132 * rSges(5,3);
t129 = t132 * pkin(6);
t213 = t155 * rSges(5,2);
t149 = t213 - t216;
t55 = -t127 - t129 + (pkin(3) - t149) * t131;
t76 = t164 * t131;
t248 = t155 * (-t179 / 0.2e1 + t174 / 0.2e1 - t177 / 0.2e1 + t172 / 0.2e1) + t156 * (-t178 / 0.2e1 - t175 / 0.2e1 - t176 / 0.2e1 - t173 / 0.2e1) - m(6) * (t161 * t50 - t51 * t76) - m(5) * (t114 * t54 - t115 * t55);
t130 = t131 ^ 2;
t246 = t132 ^ 2;
t245 = 0.2e1 * qJD(4);
t244 = 0.4e1 * qJD(4);
t243 = m(5) / 0.2e1;
t242 = m(6) / 0.2e1;
t183 = rSges(6,1) * t195 + t131 * t153 - t251;
t223 = -t129 + t251 + (pkin(3) - t258) * t131;
t159 = t132 * (t129 + (-pkin(3) - t212) * t131 + t183 + t223);
t21 = t223 * t131 + (t128 + (pkin(3) + t185) * t132 + t180) * t132;
t241 = m(6) * t159 * t21;
t191 = rSges(5,1) * t195 - t127;
t96 = t131 * t149 + t127;
t160 = t132 * (-rSges(5,2) * t196 + t191 + t96);
t39 = t132 * (-rSges(5,1) * t193 + t184) + t131 * t96;
t240 = m(5) * t39 * t160;
t236 = m(6) * (t131 * t76 + t72);
t181 = t155 * rSges(6,1) + t214;
t100 = -pkin(4) * t196 - t131 * t181;
t235 = m(6) * (t131 * t100 - t72);
t229 = m(5) * t182;
t49 = -rSges(6,2) * t196 + t183;
t227 = t49 - t51;
t53 = -t129 + (pkin(3) - t213) * t131 + t191;
t226 = t53 - t55;
t219 = m(6) * qJD(3);
t218 = m(6) * qJD(4);
t217 = m(6) * qJD(5);
t16 = t159 * t242 + t160 * t243;
t192 = t16 * qJD(4);
t189 = t257 * t233;
t188 = t257 * t234;
t187 = t169 / 0.2e1 + t171 / 0.2e1;
t186 = -t148 + t228;
t157 = (-t100 * t208 + t252) * t242 - (-t131 * t208 - t132 * t207) * t229 / 0.2e1;
t158 = (t114 * t207 + t115 * t208) * t243 + (t208 * t76 + t252) * t242;
t17 = t157 - t158;
t167 = t16 * qJD(1) + t17 * qJD(2);
t166 = -(t273 + (t260 * t156 + t269) * t131 - t266) * t131 / 0.2e1 + t265 * t234 + (t269 * t132 + t255) * t233;
t165 = (-t265 + t266 + t267) * t132 / 0.2e1 + t266 * t233 + (-t255 + t274) * t234;
t116 = -t131 * t207 + t132 * t208;
t101 = t186 * t131;
t98 = t186 * t132;
t47 = t114 * t132 + t115 * t131;
t45 = t51 * t132;
t41 = t235 / 0.2e1;
t36 = t236 / 0.2e1;
t35 = t246 * t164 + (pkin(4) * t155 + t181) * t130;
t26 = -t131 * t50 - t45;
t18 = t157 + t158;
t11 = t16 * qJD(3);
t10 = t36 - t235 / 0.2e1;
t9 = t41 + t36;
t8 = t41 - t236 / 0.2e1;
t1 = t131 * t166 + t132 * t165 + t240 + t241;
t2 = [0, 0, t192, t11 + (t35 * t242 + t47 * t243) * t245, 0; 0, 0, t18 * qJD(4) + 0.2e1 * (m(4) * ((rSges(4,1) * t131 + rSges(4,2) * t132) * t207 + (-rSges(4,1) * t132 + rSges(4,2) * t131) * t208) / 0.2e1 + (t207 * t53 + t208 * t54) * t243 + (t207 * t49 + t208 * t50) * t242) * qJD(3) + t116 * t217, t18 * qJD(3) + ((t101 * t207 - t208 * t98) * t242 + t116 * t149 * t243) * t245, t116 * t219; -t192, -t17 * qJD(4), (m(5) * t226 * t54 + m(6) * t227 * t50) * qJD(3) - t248 * qJD(4) + t26 * t217, -t248 * qJD(3) + (t101 * t50 + t51 * t98 + (-t100 - t76) * t161) * t218 + t9 * qJD(5) + (-t240 / 0.4e1 - t241 / 0.4e1) * t244 + (m(5) * (-t115 * t182 - t149 * t55) + t187 * t132 - t165) * t190 - t167 + ((m(5) * (t114 * t182 - t149 * t54) + t187 * t131 - t166) * t131 + ((t249 * t132 + t270) * t234 + (t249 * t131 - t262) * t233) * t155 + ((t250 * t132 + t260) * t234 + (t250 * t131 + t261) * t233) * t156) * qJD(4), t9 * qJD(4) + t219 * t26; t11, t17 * qJD(3), -t100 * t227 * t219 + t1 * qJD(4) + t8 * qJD(5) + t167 + (t226 * t229 * t131 + t248) * qJD(3), t1 * qJD(3) + t188 * t246 * t190 + (m(5) * (t39 * t47 - (t130 + t246) * t149 * t182) / 0.4e1 + m(6) * (-t100 * t101 + t161 * t98 + t21 * t35) / 0.4e1) * t244 + (t189 * qJD(4) * t130 + (t131 * t188 + t132 * t189) * t190) * t131, t8 * qJD(3); 0, 0, (t132 * t49 - t26 - t45) * t219 + t10 * qJD(4), t10 * qJD(3) + (t101 * t132 - t131 * t98) * t218, 0;];
Cq = t2;
