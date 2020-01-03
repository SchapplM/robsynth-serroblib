% Calculate joint inertia matrix for
% S5RRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPP1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP1_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP1_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP1_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:49:24
% EndTime: 2019-12-31 20:49:27
% DurationCPUTime: 1.04s
% Computational Cost: add. (1969->148), mult. (1604->212), div. (0->0), fcn. (1388->8), ass. (0->82)
t106 = qJ(3) + pkin(8);
t102 = cos(t106);
t180 = t102 ^ 2;
t179 = Icges(6,4) + Icges(5,5);
t178 = Icges(5,6) - Icges(6,6);
t111 = cos(qJ(3));
t174 = Icges(4,6) * t111;
t109 = sin(qJ(3));
t175 = Icges(4,5) * t109;
t101 = sin(t106);
t176 = t179 * t101;
t177 = t178 * t102;
t168 = t174 + t175 + t176 + t177;
t107 = qJ(1) + qJ(2);
t103 = sin(t107);
t99 = t103 ^ 2;
t104 = cos(t107);
t100 = t104 ^ 2;
t171 = Icges(6,2) + Icges(4,3) + Icges(5,3);
t170 = Icges(4,5) * t111 - Icges(4,6) * t109 - t178 * t101 + t179 * t102;
t165 = rSges(6,1) + pkin(4);
t164 = rSges(6,3) + qJ(5);
t166 = t164 * t101 + t165 * t102;
t163 = t171 * t103 + t170 * t104;
t162 = -t170 * t103 + t171 * t104;
t82 = t109 * rSges(4,1) + t111 * rSges(4,2);
t161 = m(4) * t82;
t158 = m(6) * t101;
t157 = pkin(3) * t109;
t110 = sin(qJ(1));
t156 = t110 * pkin(1);
t108 = -qJ(4) - pkin(7);
t97 = t111 * pkin(3) + pkin(2);
t137 = -t103 * t108 + t104 * t97;
t139 = t104 * t108;
t152 = -t104 * pkin(2) - t103 * pkin(7);
t95 = t104 * pkin(7);
t155 = t103 * (t139 + t95 + (-pkin(2) + t97) * t103) + t104 * (t137 + t152);
t154 = t103 * t101 * rSges(5,2) + t104 * rSges(5,3);
t149 = rSges(4,2) * t109;
t153 = t104 * rSges(4,3) + t103 * t149;
t151 = rSges(4,1) * t111;
t150 = rSges(5,1) * t102;
t148 = t99 + t100;
t141 = t101 * t104;
t140 = t102 * t104;
t138 = -t101 * rSges(5,1) - t102 * rSges(5,2) - t157;
t67 = t104 * rSges(3,1) - t103 * rSges(3,2);
t136 = -t165 * t101 + t164 * t102 - t157;
t135 = t103 * rSges(6,2) + t165 * t140 + t164 * t141;
t66 = -t103 * rSges(3,1) - t104 * rSges(3,2);
t116 = t103 * rSges(4,3) + (-t149 + t151) * t104;
t115 = rSges(5,1) * t140 - rSges(5,2) * t141 + t103 * rSges(5,3);
t31 = t116 - t152;
t114 = Icges(4,2) * t111 ^ 2 + Icges(3,3) + (Icges(4,1) * t109 + 0.2e1 * Icges(4,4) * t111) * t109 + (Icges(5,2) + Icges(6,3)) * t180 + ((Icges(5,1) + Icges(6,1)) * t101 + (2 * Icges(5,4) - 2 * Icges(6,5)) * t102) * t101;
t14 = t135 + t137;
t30 = t95 + (-pkin(2) - t151) * t103 + t153;
t23 = t115 + t137;
t113 = t168 * t100 + (t175 / 0.2e1 + t174 / 0.2e1 + t177 / 0.2e1 + t176 / 0.2e1 + t168 / 0.2e1) * t99;
t22 = -t139 + (-t97 - t150) * t103 + t154;
t92 = t104 * rSges(6,2);
t13 = -t139 + t92 + (-t166 - t97) * t103;
t112 = cos(qJ(1));
t105 = t112 * pkin(1);
t84 = t112 * rSges(2,1) - t110 * rSges(2,2);
t83 = -t110 * rSges(2,1) - t112 * rSges(2,2);
t56 = t105 + t67;
t55 = t66 - t156;
t45 = t138 * t104;
t44 = t138 * t103;
t27 = t105 + t31;
t26 = t30 - t156;
t25 = t136 * t104;
t24 = t136 * t103;
t21 = t105 + t23;
t20 = t22 - t156;
t17 = t103 * (t103 * t151 - t153) + t104 * t116;
t12 = t105 + t14;
t11 = t13 - t156;
t2 = t103 * (t103 * t150 - t154) + t104 * t115 + t155;
t1 = t135 * t104 + (t166 * t103 - t92) * t103 + t155;
t3 = [Icges(2,3) + m(6) * (t11 ^ 2 + t12 ^ 2) + m(5) * (t20 ^ 2 + t21 ^ 2) + m(4) * (t26 ^ 2 + t27 ^ 2) + m(3) * (t55 ^ 2 + t56 ^ 2) + m(2) * (t83 ^ 2 + t84 ^ 2) + t114; m(6) * (t13 * t11 + t14 * t12) + m(5) * (t22 * t20 + t23 * t21) + m(4) * (t30 * t26 + t31 * t27) + m(3) * (t66 * t55 + t67 * t56) + t114; m(6) * (t13 ^ 2 + t14 ^ 2) + m(5) * (t22 ^ 2 + t23 ^ 2) + m(4) * (t30 ^ 2 + t31 ^ 2) + m(3) * (t66 ^ 2 + t67 ^ 2) + t114; t113 + m(6) * (t25 * t11 + t24 * t12) + m(5) * (t45 * t20 + t44 * t21) + (-t103 * t27 - t104 * t26) * t161; t113 + m(6) * (t25 * t13 + t24 * t14) + m(5) * (t45 * t22 + t44 * t23) + (-t103 * t31 - t104 * t30) * t161; m(5) * (t2 ^ 2 + t44 ^ 2 + t45 ^ 2) + m(6) * (t1 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(4) * (t148 * t82 ^ 2 + t17 ^ 2) + t163 * t99 * t103 + (t162 * t100 + (t162 * t103 + t163 * t104) * t103) * t104; m(6) * (t103 * t11 - t104 * t12) + m(5) * (t103 * t20 - t104 * t21); m(6) * (t103 * t13 - t104 * t14) + m(5) * (t103 * t22 - t104 * t23); m(5) * (t103 * t45 - t104 * t44) + m(6) * (t103 * t25 - t104 * t24); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * t148; (t103 * t12 + t104 * t11) * t158; (t103 * t14 + t104 * t13) * t158; m(6) * (-t102 * t1 + (t103 * t24 + t104 * t25) * t101); 0; m(6) * (t148 * t101 ^ 2 + t180);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t3(1), t3(2), t3(4), t3(7), t3(11); t3(2), t3(3), t3(5), t3(8), t3(12); t3(4), t3(5), t3(6), t3(9), t3(13); t3(7), t3(8), t3(9), t3(10), t3(14); t3(11), t3(12), t3(13), t3(14), t3(15);];
Mq = res;
