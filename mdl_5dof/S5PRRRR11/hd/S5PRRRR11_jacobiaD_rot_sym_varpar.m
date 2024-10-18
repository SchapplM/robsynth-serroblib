% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRRRR11
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S5PRRRR11_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5PRRRR11_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR11_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR11_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
JaD_rot=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 21:46:36
	% EndTime: 2024-09-27 21:46:36
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	% Symbolic code from jacobiaD_rot_1_floatb_twist_matlab.m not found
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 21:46:36
	% EndTime: 2024-09-27 21:46:36
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (31->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 21:46:36
	% EndTime: 2024-09-27 21:46:37
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (573->40), mult. (853->109), div. (126->12), fcn. (1047->9), ass. (0->55)
	t102 = pkin(10) + qJ(2);
	t95 = sin(t102);
	t93 = t95 ^ 2;
	t103 = sin(pkin(5));
	t98 = t103 ^ 2;
	t133 = t93 * t98;
	t104 = cos(pkin(5));
	t96 = cos(t102);
	t124 = t96 * t103;
	t89 = atan2(t124, t104);
	t85 = sin(t89);
	t86 = cos(t89);
	t75 = t86 * t104 + t85 * t124;
	t70 = 0.1e1 / t75;
	t105 = sin(qJ(3));
	t121 = t104 * t105;
	t106 = cos(qJ(3));
	t123 = t96 * t106;
	t115 = t95 * t121 - t123;
	t78 = 0.1e1 / t115;
	t99 = 0.1e1 / t104;
	t100 = 0.1e1 / t104 ^ 2;
	t71 = 0.1e1 / t75 ^ 2;
	t79 = 0.1e1 / t115 ^ 2;
	t132 = -0.2e1 * t99 * t100;
	t82 = -t95 * t106 - t96 * t121;
	t120 = t104 * t106;
	t83 = t96 * t105 + t95 * t120;
	t69 = t82 * qJD(2) - t83 * qJD(3);
	t129 = t69 * t78 * t79;
	t117 = t96 * t120;
	t125 = t95 * t105;
	t68 = -qJD(2) * t117 - qJD(3) * t123 + (qJD(3) * t104 + qJD(2)) * t125;
	t130 = t68 * t79;
	t77 = t83 ^ 2;
	t76 = t77 * t79 + 0.1e1;
	t131 = (t77 * t129 - t83 * t130) / t76 ^ 2;
	t128 = t71 * t95;
	t127 = t82 * t83;
	t126 = t98 * t99;
	t122 = qJD(2) * t96;
	t94 = t96 ^ 2;
	t90 = t94 * t98 * t100 + 0.1e1;
	t87 = 0.1e1 / t90;
	t119 = t87 * t126;
	t88 = 0.1e1 / t90 ^ 2;
	t118 = t88 * t103 * t133;
	t116 = (t87 - 0.1e1) * t103;
	t81 = t117 - t125;
	t64 = (-t86 * t96 * t119 + t85 * t116) * t95;
	t73 = 0.1e1 / t76;
	t72 = t70 * t71;
	t67 = t71 * t133 + 0.1e1;
	t63 = qJD(2) * t64;
	t1 = [0, (-t103 * t87 * t99 + t118 * t132) * t122, 0, 0, 0; 0, (0.2e1 * (t64 * t128 - t70 * t96) / t67 ^ 2 * (-t63 * t72 * t93 + t122 * t128) * t98 + ((0.2e1 * t64 * t72 * t95 - t71 * t96) * t63 + (-t95 * t70 + ((-t64 + (-t100 * t118 - t116) * t95 * t85) * t96 - (-t94 * t119 + (0.2e1 * t119 + (t94 * t98 ^ 2 * t132 - t126) * t88) * t93) * t95 * t86) * t71) * qJD(2)) / t67) * t103, 0, 0, 0; 0, 0.2e1 * (t79 * t127 + t78 * t81) * t131 + (-(-t83 * qJD(2) + t82 * qJD(3)) * t78 - 0.2e1 * t127 * t129 + (-t81 * t69 - (t115 * qJD(2) - t81 * qJD(3)) * t83 + t82 * t68) * t79) * t73, -0.2e1 * t131 + 0.2e1 * (-t73 * t130 + (t73 * t129 - t79 * t131) * t83) * t83, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 21:46:36
	% EndTime: 2024-09-27 21:46:37
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (1878->51), mult. (1347->122), div. (158->12), fcn. (1387->13), ass. (0->65)
	t159 = pkin(10) + qJ(2);
	t151 = cos(t159);
	t161 = qJ(3) + qJ(4);
	t152 = sin(t161);
	t174 = pkin(5) + t161;
	t175 = pkin(5) - t161;
	t137 = cos(t175) / 0.2e1 + cos(t174) / 0.2e1;
	t150 = sin(t159);
	t187 = t150 * t137;
	t128 = t151 * t152 + t187;
	t122 = t128 ^ 2;
	t136 = sin(t174) / 0.2e1 - sin(t175) / 0.2e1;
	t153 = cos(t161);
	t173 = -t150 * t136 + t151 * t153;
	t124 = 0.1e1 / t173 ^ 2;
	t197 = t122 * t124;
	t148 = t150 ^ 2;
	t162 = sin(pkin(5));
	t155 = t162 ^ 2;
	t196 = t148 * t155;
	t160 = qJD(3) + qJD(4);
	t195 = qJD(2) * t152 + t136 * t160;
	t163 = cos(pkin(5));
	t186 = t151 * t162;
	t141 = atan2(t186, t163);
	t134 = sin(t141);
	t135 = cos(t141);
	t121 = t134 * t186 + t135 * t163;
	t118 = 0.1e1 / t121;
	t123 = 0.1e1 / t173;
	t156 = 0.1e1 / t163;
	t119 = 0.1e1 / t121 ^ 2;
	t157 = 0.1e1 / t163 ^ 2;
	t194 = -0.2e1 * t156 * t157;
	t114 = 0.1e1 + t197;
	t182 = qJD(2) * t151;
	t183 = t160 * t153;
	t110 = -t137 * t182 + t195 * t150 - t151 * t183;
	t190 = t124 * t128;
	t179 = t110 * t190;
	t127 = -t151 * t136 - t150 * t153;
	t132 = t137 * t160;
	t185 = t152 * t160;
	t111 = t127 * qJD(2) - t150 * t132 - t151 * t185;
	t125 = t123 * t124;
	t189 = t125 * t111;
	t193 = (-t122 * t189 - t179) / t114 ^ 2;
	t112 = 0.1e1 / t114;
	t192 = t112 * t124;
	t191 = t119 * t150;
	t184 = t155 * t156;
	t180 = 0.2e1 * t127 * t128;
	t149 = t151 ^ 2;
	t142 = t149 * t155 * t157 + 0.1e1;
	t139 = 0.1e1 / t142;
	t178 = t139 * t184;
	t140 = 0.1e1 / t142 ^ 2;
	t177 = t140 * t162 * t196;
	t176 = (t139 - 0.1e1) * t162;
	t109 = (-t135 * t151 * t178 + t134 * t176) * t150;
	t120 = t118 * t119;
	t117 = t119 * t196 + 0.1e1;
	t108 = qJD(2) * t109;
	t105 = 0.2e1 * (-t123 * t173 - t197) * t193 + (-0.2e1 * t179 + (-0.2e1 * t122 * t125 - t124 * t173 + t123) * t111) * t112;
	t1 = [0, (-t139 * t156 * t162 + t177 * t194) * t182, 0, 0, 0; 0, (0.2e1 * (t109 * t191 - t118 * t151) / t117 ^ 2 * (-t108 * t120 * t148 + t182 * t191) * t155 + ((0.2e1 * t109 * t120 * t150 - t119 * t151) * t108 + (-t150 * t118 + ((-t109 + (-t157 * t177 - t176) * t150 * t134) * t151 - (-t149 * t178 + (0.2e1 * t178 + (t149 * t155 ^ 2 * t194 - t184) * t140) * t148) * t150 * t135) * t119) * qJD(2)) / t117) * t162, 0, 0, 0; 0, t127 * t110 * t192 + t124 * t180 * t193 + (-t111 * t192 - 0.2e1 * t123 * t193) * (t151 * t137 - t150 * t152) + ((-qJD(2) * t187 - t150 * t183 - t195 * t151) * t123 - (-t173 * qJD(2) - t151 * t132 + t150 * t185) * t190 + t180 * t189) * t112, t105, t105, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 21:46:37
	% EndTime: 2024-09-27 21:46:37
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (3192->52), mult. (1659->123), div. (183->12), fcn. (1662->13), ass. (0->66)
	t166 = qJ(3) + qJ(4) + qJ(5);
	t184 = pkin(5) + t166;
	t185 = pkin(5) - t166;
	t146 = cos(t185) / 0.2e1 + cos(t184) / 0.2e1;
	t172 = pkin(10) + qJ(2);
	t163 = sin(t172);
	t159 = sin(t166);
	t164 = cos(t172);
	t197 = t164 * t159;
	t139 = t163 * t146 + t197;
	t133 = t139 ^ 2;
	t145 = sin(t184) / 0.2e1 - sin(t185) / 0.2e1;
	t160 = cos(t166);
	t186 = -t163 * t145 + t164 * t160;
	t135 = 0.1e1 / t186 ^ 2;
	t207 = t133 * t135;
	t161 = t163 ^ 2;
	t173 = sin(pkin(5));
	t168 = t173 ^ 2;
	t206 = t161 * t168;
	t174 = cos(pkin(5));
	t196 = t164 * t173;
	t152 = atan2(t196, t174);
	t148 = sin(t152);
	t149 = cos(t152);
	t132 = t148 * t196 + t149 * t174;
	t129 = 0.1e1 / t132;
	t134 = 0.1e1 / t186;
	t169 = 0.1e1 / t174;
	t130 = 0.1e1 / t132 ^ 2;
	t170 = 0.1e1 / t174 ^ 2;
	t205 = -0.2e1 * t169 * t170;
	t125 = 0.1e1 + t207;
	t165 = qJD(3) + qJD(4) + qJD(5);
	t183 = t145 * t165;
	t192 = qJD(2) * t164;
	t193 = qJD(2) * t163;
	t195 = t165 * t160;
	t121 = -t146 * t192 + t159 * t193 + t163 * t183 - t164 * t195;
	t201 = t135 * t139;
	t190 = t121 * t201;
	t138 = -t164 * t145 - t163 * t160;
	t143 = t146 * t165;
	t122 = t138 * qJD(2) - t163 * t143 - t165 * t197;
	t136 = t134 * t135;
	t200 = t136 * t122;
	t204 = (-t133 * t200 - t190) / t125 ^ 2;
	t123 = 0.1e1 / t125;
	t203 = t123 * t135;
	t202 = t130 * t163;
	t198 = t163 * t159;
	t194 = t168 * t169;
	t191 = 0.2e1 * t138 * t139;
	t162 = t164 ^ 2;
	t153 = t162 * t168 * t170 + 0.1e1;
	t150 = 0.1e1 / t153;
	t189 = t150 * t194;
	t151 = 0.1e1 / t153 ^ 2;
	t188 = t151 * t173 * t206;
	t187 = (t150 - 0.1e1) * t173;
	t120 = (-t149 * t164 * t189 + t148 * t187) * t163;
	t131 = t129 * t130;
	t128 = t130 * t206 + 0.1e1;
	t119 = qJD(2) * t120;
	t116 = 0.2e1 * (-t134 * t186 - t207) * t204 + (-0.2e1 * t190 + (-0.2e1 * t133 * t136 - t135 * t186 + t134) * t122) * t123;
	t1 = [0, (-t150 * t169 * t173 + t188 * t205) * t192, 0, 0, 0; 0, (0.2e1 * (t120 * t202 - t129 * t164) / t128 ^ 2 * (-t119 * t131 * t161 + t192 * t202) * t168 + ((0.2e1 * t120 * t131 * t163 - t130 * t164) * t119 + (-t163 * t129 + ((-t120 + (-t170 * t188 - t187) * t163 * t148) * t164 - (-t162 * t189 + (0.2e1 * t189 + (t162 * t168 ^ 2 * t205 - t194) * t151) * t161) * t163 * t149) * t130) * qJD(2)) / t128) * t173, 0, 0, 0; 0, t138 * t121 * t203 + t135 * t191 * t204 + (-t122 * t203 - 0.2e1 * t134 * t204) * (t164 * t146 - t198) + ((-t146 * t193 - t159 * t192 - t163 * t195 - t164 * t183) * t134 - (-t186 * qJD(2) - t164 * t143 + t165 * t198) * t201 + t191 * t200) * t123, t116, t116, t116;];
	JaD_rot = t1;
end