% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRPPRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPPRR1_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR1_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_jacobiaD_rot_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:44:35
% EndTime: 2019-02-26 19:44:35
% DurationCPUTime: 0.40s
% Computational Cost: add. (1475->50), mult. (4551->125), div. (254->12), fcn. (5909->13), ass. (0->70)
t150 = cos(pkin(6));
t144 = sin(pkin(11));
t148 = cos(pkin(11));
t151 = sin(qJ(2));
t152 = cos(qJ(2));
t162 = t152 * t144 + t151 * t148;
t136 = t162 * t150;
t139 = t151 * t144 - t152 * t148;
t137 = t139 * qJD(2);
t145 = sin(pkin(10));
t149 = cos(pkin(10));
t160 = t139 * t150;
t121 = -t145 * t162 - t149 * t160;
t146 = sin(pkin(6));
t134 = t139 * t146;
t107 = atan2(t121, t134);
t102 = sin(t107);
t103 = cos(t107);
t98 = t102 * t121 + t103 * t134;
t95 = 0.1e1 / t98;
t143 = sin(pkin(12));
t147 = cos(pkin(12));
t163 = -t145 * t136 - t149 * t139;
t168 = t145 * t146;
t113 = t143 * t168 + t147 * t163;
t109 = 0.1e1 / t113;
t131 = 0.1e1 / t134;
t110 = 0.1e1 / t113 ^ 2;
t132 = 0.1e1 / t134 ^ 2;
t96 = 0.1e1 / t98 ^ 2;
t112 = t143 * t163 - t147 * t168;
t108 = t112 ^ 2;
t101 = t108 * t110 + 0.1e1;
t111 = t109 * t110;
t130 = t150 * t137;
t138 = t162 * qJD(2);
t117 = t145 * t130 - t149 * t138;
t176 = 0.1e1 / t101 ^ 2 * (-t108 * t111 * t147 + t110 * t112 * t143) * t117;
t99 = 0.1e1 / t101;
t175 = t117 * t99;
t123 = t145 * t160 - t149 * t162;
t174 = t123 * t96;
t173 = t109 * t143;
t172 = t112 * t147;
t171 = t121 * t132;
t135 = t162 * t146;
t170 = t121 * t135;
t128 = qJD(2) * t135;
t169 = t128 * t131 * t132;
t120 = -t149 * t136 + t145 * t139;
t118 = t121 ^ 2;
t106 = t118 * t132 + 0.1e1;
t104 = 0.1e1 / t106;
t161 = -t120 * t131 + t132 * t170;
t90 = t161 * t104;
t165 = t134 * t90 + t120;
t164 = -t121 * t90 + t135;
t159 = qJD(2) * t136;
t129 = t146 * t137;
t119 = t123 ^ 2;
t116 = t149 * t137 + t145 * t159;
t115 = t149 * t130 + t145 * t138;
t114 = t145 * t137 - t149 * t159;
t97 = t95 * t96;
t93 = t119 * t96 + 0.1e1;
t89 = (t114 * t131 - t128 * t171) * t104;
t87 = t165 * t102 + t164 * t103;
t86 = (t121 * t89 + t128) * t103 + (-t134 * t89 + t114) * t102;
t85 = 0.2e1 * t161 * (t114 * t171 - t118 * t169) / t106 ^ 2 + (0.2e1 * t169 * t170 + t115 * t131 + (-t114 * t135 - t120 * t128 + t121 * t129) * t132) * t104;
t1 = [0, t85, 0, 0, 0, 0; 0, 0.2e1 * (-t163 * t95 - t87 * t174) / t93 ^ 2 * (-t119 * t97 * t86 + t116 * t174) + (t117 * t95 + t87 * t116 * t96 + (-0.2e1 * t87 * t123 * t97 - t163 * t96) * t86 + ((-t114 * t90 + t121 * t85 + t165 * t89 - t129) * t103 + (t128 * t90 - t134 * t85 - t164 * t89 + t115) * t102) * t174) / t93, 0, 0, 0, 0; 0 (-t110 * t172 + t173) * t99 * t116 + 0.2e1 * (-t173 * t176 + (t111 * t172 * t175 + (t112 * t176 - t143 * t175) * t110) * t147) * t123, 0, 0, 0, 0;];
JaD_rot  = t1;
