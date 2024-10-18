% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRRRR14
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
%   Wie in S5RRRRR14_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RRRRR14_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR14_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR14_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
JaD_rot=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 18:44:41
	% EndTime: 2024-09-27 18:44:41
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 18:44:41
	% EndTime: 2024-09-27 18:44:41
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 18:44:42
	% EndTime: 2024-09-27 18:44:42
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (74->0), mult. (74->0), div. (30->0), fcn. (44->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 18:44:42
	% EndTime: 2024-09-27 18:44:42
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (1102->41), mult. (1512->108), div. (234->12), fcn. (1884->9), ass. (0->60)
	t123 = qJ(1) + qJ(2);
	t115 = sin(t123);
	t113 = t115 ^ 2;
	t124 = sin(pkin(5));
	t118 = t124 ^ 2;
	t155 = t113 * t118;
	t125 = cos(pkin(5));
	t116 = cos(t123);
	t145 = t116 * t124;
	t110 = atan2(t145, t125);
	t105 = sin(t110);
	t106 = cos(t110);
	t95 = t105 * t145 + t106 * t125;
	t90 = 0.1e1 / t95;
	t126 = sin(qJ(3));
	t142 = t125 * t126;
	t127 = cos(qJ(3));
	t144 = t116 * t127;
	t136 = t115 * t142 - t144;
	t98 = 0.1e1 / t136;
	t119 = 0.1e1 / t125;
	t99 = 0.1e1 / t136 ^ 2;
	t120 = 0.1e1 / t125 ^ 2;
	t91 = 0.1e1 / t95 ^ 2;
	t154 = -0.2e1 * t119 * t120;
	t141 = t125 * t127;
	t103 = t115 * t141 + t116 * t126;
	t102 = -t115 * t127 - t116 * t142;
	t122 = qJD(1) + qJD(2);
	t89 = -t103 * qJD(3) + t102 * t122;
	t151 = t98 * t99 * t89;
	t138 = t116 * t141;
	t147 = t115 * t126;
	t88 = -t122 * t138 - qJD(3) * t144 + (qJD(3) * t125 + t122) * t147;
	t152 = t88 * t99;
	t97 = t103 ^ 2;
	t96 = t97 * t99 + 0.1e1;
	t153 = (-t103 * t152 + t97 * t151) / t96 ^ 2;
	t150 = t115 * t91;
	t149 = t116 * t91;
	t148 = t102 * t103;
	t146 = t116 * t122;
	t143 = t118 * t119;
	t114 = t116 ^ 2;
	t111 = t114 * t118 * t120 + 0.1e1;
	t108 = 0.1e1 / t111;
	t140 = t108 * t143;
	t109 = 0.1e1 / t111 ^ 2;
	t139 = t109 * t124 * t155;
	t137 = (t108 - 0.1e1) * t124;
	t101 = t138 - t147;
	t83 = (-t106 * t116 * t140 + t105 * t137) * t115;
	t93 = 0.1e1 / t96;
	t92 = t90 * t91;
	t87 = t91 * t155 + 0.1e1;
	t84 = (-t108 * t119 * t124 + t139 * t154) * t146;
	t82 = t122 * t83;
	t79 = 0.2e1 * (t101 * t98 + t99 * t148) * t153 + (-(t102 * qJD(3) - t103 * t122) * t98 - 0.2e1 * t148 * t151 + (-t101 * t89 - (-t101 * qJD(3) + t136 * t122) * t103 + t102 * t88) * t99) * t93;
	t78 = (0.2e1 * (-t116 * t90 + t83 * t150) / t87 ^ 2 * (-t113 * t82 * t92 + t146 * t150) * t118 + ((0.2e1 * t115 * t83 * t92 - t149) * t82 + (-t83 * t149 + (-t90 + (-t120 * t139 - t137) * t105 * t149 - (-t114 * t140 + (0.2e1 * t140 + (t114 * t118 ^ 2 * t154 - t143) * t109) * t113) * t91 * t106) * t115) * t122) / t87) * t124;
	t1 = [t84, t84, 0, 0, 0; t78, t78, 0, 0, 0; t79, t79, -0.2e1 * t153 + 0.2e1 * (-t93 * t152 + (t93 * t151 - t99 * t153) * t103) * t103, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 18:44:42
	% EndTime: 2024-09-27 18:44:42
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (2786->53), mult. (2070->124), div. (266->12), fcn. (2224->13), ass. (0->69)
	t184 = qJ(3) + qJ(4);
	t200 = pkin(5) + t184;
	t201 = pkin(5) - t184;
	t160 = cos(t201) / 0.2e1 + cos(t200) / 0.2e1;
	t173 = sin(t184);
	t185 = qJ(1) + qJ(2);
	t174 = sin(t185);
	t176 = cos(t185);
	t151 = t160 * t174 + t173 * t176;
	t145 = t151 ^ 2;
	t159 = sin(t200) / 0.2e1 - sin(t201) / 0.2e1;
	t175 = cos(t184);
	t198 = -t159 * t174 + t176 * t175;
	t147 = 0.1e1 / t198 ^ 2;
	t219 = t145 * t147;
	t171 = t174 ^ 2;
	t186 = sin(pkin(5));
	t178 = t186 ^ 2;
	t218 = t171 * t178;
	t187 = cos(pkin(5));
	t210 = t176 * t186;
	t164 = atan2(t210, t187);
	t157 = sin(t164);
	t158 = cos(t164);
	t144 = t157 * t210 + t158 * t187;
	t141 = 0.1e1 / t144;
	t146 = 0.1e1 / t198;
	t179 = 0.1e1 / t187;
	t142 = 0.1e1 / t144 ^ 2;
	t180 = 0.1e1 / t187 ^ 2;
	t217 = -0.2e1 * t179 * t180;
	t136 = 0.1e1 + t219;
	t182 = qJD(3) + qJD(4);
	t196 = t159 * t182;
	t183 = qJD(1) + qJD(2);
	t207 = t183 * t176;
	t208 = t182 * t175;
	t211 = t174 * t183;
	t132 = -t160 * t207 + t173 * t211 + t174 * t196 - t176 * t208;
	t213 = t147 * t151;
	t205 = t132 * t213;
	t197 = t159 * t183 + t173 * t182;
	t199 = -t160 * t182 - t175 * t183;
	t133 = t174 * t199 - t176 * t197;
	t148 = t146 * t147;
	t212 = t148 * t133;
	t216 = (-t145 * t212 - t205) / t136 ^ 2;
	t134 = 0.1e1 / t136;
	t215 = t134 * t147;
	t214 = t142 * t174;
	t209 = t178 * t179;
	t150 = -t159 * t176 - t174 * t175;
	t206 = 0.2e1 * t150 * t151;
	t172 = t176 ^ 2;
	t165 = t172 * t178 * t180 + 0.1e1;
	t162 = 0.1e1 / t165;
	t204 = t162 * t209;
	t163 = 0.1e1 / t165 ^ 2;
	t203 = t163 * t186 * t218;
	t202 = (t162 - 0.1e1) * t186;
	t131 = (-t158 * t176 * t204 + t157 * t202) * t174;
	t143 = t141 * t142;
	t140 = t142 * t218 + 0.1e1;
	t137 = (-t162 * t179 * t186 + t203 * t217) * t207;
	t130 = t183 * t131;
	t127 = t150 * t132 * t215 + t147 * t206 * t216 + (-t133 * t215 - 0.2e1 * t146 * t216) * (t160 * t176 - t173 * t174) + ((-t160 * t211 - t173 * t207 - t174 * t208 - t176 * t196) * t146 - (t174 * t197 + t176 * t199) * t213 + t206 * t212) * t134;
	t126 = 0.2e1 * (-t146 * t198 - t219) * t216 + (-0.2e1 * t205 + (-0.2e1 * t145 * t148 - t147 * t198 + t146) * t133) * t134;
	t125 = (0.2e1 * (t131 * t214 - t141 * t176) / t140 ^ 2 * (-t130 * t143 * t171 + t207 * t214) * t178 + ((0.2e1 * t131 * t143 * t174 - t142 * t176) * t130 + (-t174 * t141 + ((-t131 + (-t180 * t203 - t202) * t174 * t157) * t176 - (-t172 * t204 + (0.2e1 * t204 + (t172 * t178 ^ 2 * t217 - t209) * t163) * t171) * t174 * t158) * t142) * t183) / t140) * t186;
	t1 = [t137, t137, 0, 0, 0; t125, t125, 0, 0, 0; t127, t127, t126, t126, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 18:44:42
	% EndTime: 2024-09-27 18:44:42
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (4310->53), mult. (2382->123), div. (291->12), fcn. (2499->13), ass. (0->70)
	t188 = qJ(3) + qJ(4) + qJ(5);
	t208 = pkin(5) + t188;
	t209 = pkin(5) - t188;
	t168 = cos(t209) / 0.2e1 + cos(t208) / 0.2e1;
	t181 = sin(t188);
	t195 = qJ(1) + qJ(2);
	t186 = sin(t195);
	t187 = cos(t195);
	t161 = t168 * t186 + t181 * t187;
	t155 = t161 ^ 2;
	t167 = sin(t208) / 0.2e1 - sin(t209) / 0.2e1;
	t182 = cos(t188);
	t176 = t187 * t182;
	t210 = -t167 * t186 + t176;
	t157 = 0.1e1 / t210 ^ 2;
	t229 = t155 * t157;
	t183 = t186 ^ 2;
	t196 = sin(pkin(5));
	t190 = t196 ^ 2;
	t228 = t183 * t190;
	t197 = cos(pkin(5));
	t219 = t187 * t196;
	t174 = atan2(t219, t197);
	t170 = sin(t174);
	t171 = cos(t174);
	t154 = t170 * t219 + t171 * t197;
	t151 = 0.1e1 / t154;
	t156 = 0.1e1 / t210;
	t191 = 0.1e1 / t197;
	t152 = 0.1e1 / t154 ^ 2;
	t192 = 0.1e1 / t197 ^ 2;
	t227 = -0.2e1 * t191 * t192;
	t146 = 0.1e1 + t229;
	t185 = qJD(3) + qJD(4) + qJD(5);
	t206 = t167 * t185;
	t194 = qJD(1) + qJD(2);
	t217 = t194 * t187;
	t220 = t186 * t194;
	t142 = -t168 * t217 - t176 * t185 + t181 * t220 + t186 * t206;
	t223 = t157 * t161;
	t215 = t142 * t223;
	t207 = t167 * t194 + t181 * t185;
	t211 = -t168 * t185 - t182 * t194;
	t143 = t186 * t211 - t187 * t207;
	t158 = t156 * t157;
	t222 = t158 * t143;
	t226 = (-t155 * t222 - t215) / t146 ^ 2;
	t144 = 0.1e1 / t146;
	t225 = t144 * t157;
	t224 = t152 * t186;
	t221 = t186 * t182;
	t218 = t190 * t191;
	t160 = -t167 * t187 - t221;
	t216 = 0.2e1 * t160 * t161;
	t184 = t187 ^ 2;
	t175 = t184 * t190 * t192 + 0.1e1;
	t172 = 0.1e1 / t175;
	t214 = t172 * t218;
	t173 = 0.1e1 / t175 ^ 2;
	t213 = t173 * t196 * t228;
	t212 = (t172 - 0.1e1) * t196;
	t141 = (-t171 * t187 * t214 + t170 * t212) * t186;
	t153 = t151 * t152;
	t150 = t152 * t228 + 0.1e1;
	t147 = (-t172 * t191 * t196 + t213 * t227) * t217;
	t140 = t194 * t141;
	t137 = (0.2e1 * (t141 * t224 - t151 * t187) / t150 ^ 2 * (-t140 * t153 * t183 + t217 * t224) * t190 + ((0.2e1 * t141 * t153 * t186 - t152 * t187) * t140 + (-t186 * t151 + ((-t141 + (-t192 * t213 - t212) * t186 * t170) * t187 - (-t184 * t214 + (0.2e1 * t214 + (t184 * t190 ^ 2 * t227 - t218) * t173) * t183) * t186 * t171) * t152) * t194) / t150) * t196;
	t136 = t160 * t142 * t225 + t157 * t216 * t226 + (-t143 * t225 - 0.2e1 * t156 * t226) * (t168 * t187 - t181 * t186) + ((-t168 * t220 - t181 * t217 - t185 * t221 - t187 * t206) * t156 - (t186 * t207 + t187 * t211) * t223 + t216 * t222) * t144;
	t135 = 0.2e1 * (-t156 * t210 - t229) * t226 + (-0.2e1 * t215 + (-0.2e1 * t155 * t158 - t157 * t210 + t156) * t143) * t144;
	t1 = [t147, t147, 0, 0, 0; t137, t137, 0, 0, 0; t136, t136, t135, t135, t135;];
	JaD_rot = t1;
end