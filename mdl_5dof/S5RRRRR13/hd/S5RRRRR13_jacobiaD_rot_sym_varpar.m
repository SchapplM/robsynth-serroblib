% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRRRR13
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
%   Wie in S5RRRRR13_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RRRRR13_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR13_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR13_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
JaD_rot=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 17:33:11
	% EndTime: 2024-09-27 17:33:11
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 17:33:11
	% EndTime: 2024-09-27 17:33:11
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 17:33:11
	% EndTime: 2024-09-27 17:33:11
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (74->0), mult. (74->0), div. (30->0), fcn. (44->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 17:33:11
	% EndTime: 2024-09-27 17:33:11
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (195->0), mult. (111->0), div. (45->0), fcn. (66->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 17:33:11
	% EndTime: 2024-09-27 17:33:11
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (2623->41), mult. (2171->108), div. (342->12), fcn. (2721->9), ass. (0->60)
	t123 = qJ(1) + qJ(2) + qJ(3);
	t120 = sin(t123);
	t118 = t120 ^ 2;
	t129 = sin(pkin(5));
	t125 = t129 ^ 2;
	t160 = t118 * t125;
	t130 = cos(pkin(5));
	t121 = cos(t123);
	t150 = t121 * t129;
	t115 = atan2(t150, t130);
	t111 = sin(t115);
	t112 = cos(t115);
	t101 = t111 * t150 + t112 * t130;
	t98 = 0.1e1 / t101;
	t131 = sin(qJ(4));
	t147 = t130 * t131;
	t132 = cos(qJ(4));
	t149 = t121 * t132;
	t141 = t120 * t147 - t149;
	t103 = 0.1e1 / t141;
	t126 = 0.1e1 / t130;
	t99 = 0.1e1 / t101 ^ 2;
	t104 = 0.1e1 / t141 ^ 2;
	t127 = 0.1e1 / t130 ^ 2;
	t159 = -0.2e1 * t126 * t127;
	t146 = t130 * t132;
	t108 = t120 * t146 + t121 * t131;
	t102 = t108 ^ 2;
	t122 = qJD(1) + qJD(2) + qJD(3);
	t143 = t121 * t146;
	t152 = t120 * t131;
	t93 = -t122 * t143 - qJD(4) * t149 + (qJD(4) * t130 + t122) * t152;
	t154 = t93 * t104;
	t107 = -t120 * t132 - t121 * t147;
	t94 = -qJD(4) * t108 + t107 * t122;
	t157 = t103 * t104 * t94;
	t97 = t102 * t104 + 0.1e1;
	t158 = (t102 * t157 - t108 * t154) / t97 ^ 2;
	t156 = t120 * t99;
	t155 = t121 * t99;
	t153 = t107 * t108;
	t151 = t121 * t122;
	t148 = t125 * t126;
	t119 = t121 ^ 2;
	t116 = t119 * t125 * t127 + 0.1e1;
	t113 = 0.1e1 / t116;
	t145 = t113 * t148;
	t114 = 0.1e1 / t116 ^ 2;
	t144 = t114 * t129 * t160;
	t142 = (t113 - 0.1e1) * t129;
	t106 = t143 - t152;
	t88 = (-t112 * t121 * t145 + t111 * t142) * t120;
	t100 = t98 * t99;
	t95 = 0.1e1 / t97;
	t92 = t99 * t160 + 0.1e1;
	t89 = (-t113 * t126 * t129 + t144 * t159) * t151;
	t87 = t122 * t88;
	t84 = 0.2e1 * (t103 * t106 + t104 * t153) * t158 + (-(qJD(4) * t107 - t108 * t122) * t103 - 0.2e1 * t153 * t157 + (-t106 * t94 - (-qJD(4) * t106 + t122 * t141) * t108 + t107 * t93) * t104) * t95;
	t83 = (0.2e1 * (-t121 * t98 + t156 * t88) / t92 ^ 2 * (-t100 * t118 * t87 + t151 * t156) * t125 + ((0.2e1 * t100 * t120 * t88 - t155) * t87 + (-t88 * t155 + (-t98 + (-t127 * t144 - t142) * t111 * t155 - (-t119 * t145 + (0.2e1 * t145 + (t119 * t125 ^ 2 * t159 - t148) * t114) * t118) * t99 * t112) * t120) * t122) / t92) * t129;
	t1 = [t89, t89, t89, 0, 0; t83, t83, t83, 0, 0; t84, t84, t84, -0.2e1 * t158 + 0.2e1 * (-t95 * t154 + (-t104 * t158 + t95 * t157) * t108) * t108, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 17:33:11
	% EndTime: 2024-09-27 17:33:11
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (4822->52), mult. (2793->122), div. (374->12), fcn. (3061->13), ass. (0->69)
	t183 = qJ(1) + qJ(2) + qJ(3);
	t179 = cos(t183);
	t190 = qJ(4) + qJ(5);
	t181 = sin(t190);
	t205 = pkin(5) + t190;
	t206 = pkin(5) - t190;
	t168 = cos(t206) / 0.2e1 + cos(t205) / 0.2e1;
	t178 = sin(t183);
	t217 = t168 * t178;
	t156 = t179 * t181 + t217;
	t150 = t156 ^ 2;
	t167 = sin(t205) / 0.2e1 - sin(t206) / 0.2e1;
	t182 = cos(t190);
	t203 = -t167 * t178 + t179 * t182;
	t152 = 0.1e1 / t203 ^ 2;
	t226 = t150 * t152;
	t174 = t178 ^ 2;
	t191 = sin(pkin(5));
	t185 = t191 ^ 2;
	t225 = t174 * t185;
	t180 = qJD(1) + qJD(2) + qJD(3);
	t189 = qJD(4) + qJD(5);
	t224 = t189 * t167 + t180 * t181;
	t192 = cos(pkin(5));
	t215 = t179 * t191;
	t169 = atan2(t215, t192);
	t162 = sin(t169);
	t163 = cos(t169);
	t149 = t162 * t215 + t163 * t192;
	t146 = 0.1e1 / t149;
	t151 = 0.1e1 / t203;
	t186 = 0.1e1 / t192;
	t147 = 0.1e1 / t149 ^ 2;
	t187 = 0.1e1 / t192 ^ 2;
	t223 = -0.2e1 * t186 * t187;
	t142 = 0.1e1 + t226;
	t213 = t182 * t189;
	t216 = t179 * t180;
	t137 = -t168 * t216 + t224 * t178 - t179 * t213;
	t218 = t152 * t156;
	t210 = t137 * t218;
	t202 = t167 * t180 + t181 * t189;
	t204 = -t168 * t189 - t180 * t182;
	t138 = t204 * t178 - t202 * t179;
	t153 = t151 * t152;
	t221 = t138 * t153;
	t222 = (-t150 * t221 - t210) / t142 ^ 2;
	t140 = 0.1e1 / t142;
	t220 = t140 * t152;
	t219 = t147 * t178;
	t212 = t185 * t186;
	t155 = -t167 * t179 - t178 * t182;
	t211 = 0.2e1 * t155 * t156;
	t175 = t179 ^ 2;
	t170 = t175 * t185 * t187 + 0.1e1;
	t165 = 0.1e1 / t170;
	t209 = t165 * t212;
	t166 = 0.1e1 / t170 ^ 2;
	t208 = t166 * t191 * t225;
	t207 = (t165 - 0.1e1) * t191;
	t136 = (-t163 * t179 * t209 + t162 * t207) * t178;
	t148 = t146 * t147;
	t145 = t147 * t225 + 0.1e1;
	t139 = (-t165 * t186 * t191 + t208 * t223) * t216;
	t135 = t180 * t136;
	t132 = t155 * t137 * t220 + t152 * t211 * t222 + (-t138 * t220 - 0.2e1 * t151 * t222) * (t168 * t179 - t178 * t181) + ((-t178 * t213 - t224 * t179 - t180 * t217) * t151 - (t202 * t178 + t204 * t179) * t218 + t211 * t221) * t140;
	t131 = 0.2e1 * (-t151 * t203 - t226) * t222 + (-0.2e1 * t210 + (-0.2e1 * t153 * t150 - t152 * t203 + t151) * t138) * t140;
	t130 = (0.2e1 * (t136 * t219 - t146 * t179) / t145 ^ 2 * (-t135 * t148 * t174 + t216 * t219) * t185 + ((0.2e1 * t136 * t148 * t178 - t147 * t179) * t135 + (-t178 * t146 + ((-t136 + (-t187 * t208 - t207) * t178 * t162) * t179 - (-t175 * t209 + (0.2e1 * t209 + (t175 * t185 ^ 2 * t223 - t212) * t166) * t174) * t178 * t163) * t147) * t180) / t145) * t191;
	t1 = [t139, t139, t139, 0, 0; t130, t130, t130, 0, 0; t132, t132, t132, t131, t131;];
	JaD_rot = t1;
end