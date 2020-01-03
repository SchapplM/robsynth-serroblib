% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
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
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPPR5_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR5_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR5_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR5_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR5_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:46:21
% EndTime: 2019-12-31 17:46:24
% DurationCPUTime: 1.42s
% Computational Cost: add. (3949->143), mult. (6395->226), div. (0->0), fcn. (7346->8), ass. (0->101)
t165 = sin(pkin(7));
t166 = cos(pkin(7));
t180 = sin(qJ(1));
t181 = cos(qJ(1));
t118 = t165 * t181 - t166 * t180;
t116 = t118 ^ 2;
t117 = -t165 * t180 - t166 * t181;
t129 = pkin(8) + qJ(5);
t126 = sin(t129);
t127 = cos(t129);
t141 = Icges(6,5) * t127 - Icges(6,6) * t126;
t66 = Icges(6,3) * t118 - t117 * t141;
t217 = t118 * t66;
t95 = t117 * t180 - t118 * t181;
t216 = (m(5) + m(6)) * t95;
t215 = t95 * (-m(5) / 0.2e1 - m(6) / 0.2e1);
t163 = Icges(6,4) * t127;
t143 = -Icges(6,2) * t126 + t163;
t69 = Icges(6,6) * t118 - t117 * t143;
t169 = t126 * t69;
t64 = Icges(6,3) * t117 + t118 * t141;
t214 = t64 - t169;
t146 = rSges(6,1) * t126 + rSges(6,2) * t127;
t205 = t117 ^ 2;
t211 = t146 * (t116 + t205);
t130 = cos(pkin(8));
t125 = pkin(4) * t130 + pkin(3);
t148 = -pkin(1) * t180 + qJ(2) * t181;
t136 = -pkin(2) * t180 + t148;
t171 = rSges(6,1) * t127;
t147 = -rSges(6,2) * t126 + t171;
t167 = rSges(6,3) + pkin(6) + qJ(4);
t133 = t167 * t117 + (t125 + t147) * t118 + t136;
t164 = Icges(6,4) * t126;
t142 = Icges(6,2) * t127 + t164;
t145 = Icges(6,1) * t127 - t164;
t158 = t117 * t126;
t100 = rSges(6,2) * t158;
t137 = pkin(1) * t181 + qJ(2) * t180;
t135 = pkin(2) * t181 + t137;
t41 = t100 + t167 * t118 + (-t125 - t171) * t117 + t135;
t87 = t146 * t118;
t88 = t146 * t117;
t210 = (t145 / 0.2e1 - t142 / 0.2e1) * t126 + m(6) * (-t133 * t87 + t41 * t88);
t208 = rSges(5,3) + qJ(4);
t207 = sin(pkin(8)) * rSges(5,2) - rSges(5,1) * t130 - pkin(3);
t67 = Icges(6,6) * t117 + t118 * t143;
t70 = Icges(6,5) * t117 + t118 * t145;
t203 = 4 * qJD(1);
t195 = m(3) * ((rSges(3,3) * t181 + t148) * t181 + (rSges(3,3) * t180 + t137) * t180);
t194 = m(4) * ((rSges(4,1) * t118 - rSges(4,2) * t117 + t136) * t181 + (-rSges(4,1) * t117 - rSges(4,2) * t118 + t135) * t180);
t132 = t117 * t208 - t118 * t207 + t136;
t47 = t117 * t207 + t118 * t208 + t135;
t193 = m(5) * (t117 * t132 + t118 * t47);
t192 = m(5) * (t132 * t181 + t180 * t47);
t11 = t117 * t133 + t118 * t41;
t190 = m(6) * t11;
t189 = m(6) * (t133 * t181 + t180 * t41);
t188 = m(6) * (-t180 * t87 - t181 * t88);
t186 = m(6) * (t117 * t181 + t118 * t180) * t146;
t185 = -t117 / 0.2e1;
t183 = t118 / 0.2e1;
t72 = Icges(6,5) * t118 - t117 * t145;
t168 = t127 * t72;
t177 = -t117 * t66 - t118 * t168;
t157 = t117 * t127;
t176 = -t157 * t72 + t217;
t175 = t216 / 0.2e1;
t173 = -t216 / 0.2e1;
t172 = m(6) * qJD(5);
t159 = t117 * t118;
t35 = t117 * t88 + t118 * t87;
t138 = m(6) * t35;
t139 = m(6) * t211 / 0.2e1;
t25 = t138 / 0.2e1 + t139;
t156 = t25 * qJD(1);
t144 = Icges(6,1) * t126 + t163;
t140 = Icges(6,5) * t126 + Icges(6,6) * t127;
t20 = -t118 * t64 + t157 * t70 - t158 * t67;
t82 = t140 * t117;
t81 = t140 * t118;
t61 = t186 / 0.2e1;
t50 = t188 / 0.2e1;
t26 = -t138 / 0.2e1 + t139;
t21 = t158 * t69 + t176;
t19 = t118 * t169 + t177;
t17 = t50 + t61;
t16 = t50 - t186 / 0.2e1;
t15 = t61 - t188 / 0.2e1;
t14 = t173 + t215;
t13 = t173 + t175;
t12 = t175 - t215;
t8 = (t144 / 0.2e1 + t143 / 0.2e1) * t127 + t210;
t7 = t190 + t193;
t6 = -t117 * t20 + t118 * t21;
t5 = -t117 * (-(t126 * t67 - t127 * t70) * t118 + t117 * t64) + t118 * t19;
t4 = t189 + t192 + t194 + t195;
t3 = (t19 - t20 - t177) * t117 + t176 * t118;
t2 = (t20 + (t168 + t214) * t118) * t118 + (t117 * t214 - t176 + t21 + t217) * t117;
t1 = (t2 / 0.2e1 + t5 / 0.2e1) * t118 + (t6 / 0.2e1 - t3 / 0.2e1) * t117;
t9 = [t4 * qJD(2) + t7 * qJD(4) + t8 * qJD(5), qJD(1) * t4 + qJD(4) * t13 + qJD(5) * t17, 0, t7 * qJD(1) + t13 * qJD(2) + t26 * qJD(5), t8 * qJD(1) + t17 * qJD(2) + t26 * qJD(4) + (((-t117 * t142 - t72) * t127 + (-t117 * t144 + t69) * t126) * t183 + t117 * t3 / 0.2e1 + (t11 * t147 - (t117 * t87 - t118 * t88) * t146) * m(6) - (t116 / 0.2e1 + t205 / 0.2e1) * t141 + ((-t118 * t142 + t70) * t127 + (-t118 * t144 - t67) * t126 + t6) * t185 - (t2 + t5) * t118 / 0.2e1) * qJD(5); t12 * qJD(4) + t16 * qJD(5) + (-t195 / 0.4e1 - t194 / 0.4e1 - t192 / 0.4e1 - t189 / 0.4e1) * t203, 0, 0, t12 * qJD(1), t147 * t172 * t95 + t16 * qJD(1); 0, 0, 0, 0, -t35 * t172; t14 * qJD(2) - t25 * qJD(5) + (-t193 / 0.4e1 - t190 / 0.4e1) * t203, t14 * qJD(1), 0, 0, -t156; t15 * qJD(2) + t25 * qJD(4) + t1 * qJD(5) + ((-t143 - t144) * t127 / 0.2e1 - t210) * qJD(1), t15 * qJD(1), 0, t156, t1 * qJD(1) + (m(6) * (t147 * t211 + (-t116 * t147 + t117 * (-rSges(6,1) * t157 + t100)) * t35) + (t116 * t82 - t159 * t81) * t183 + (-t159 * t82 + t205 * t81) * t185) * qJD(5);];
Cq = t9;
