% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPPRRR4
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
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S6RPPRRR4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:06
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRRR4_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR4_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRR4_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:06:31
	% EndTime: 2019-10-10 00:06:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:06:31
	% EndTime: 2019-10-10 00:06:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:06:31
	% EndTime: 2019-10-10 00:06:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:06:31
	% EndTime: 2019-10-10 00:06:32
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (52->9), mult. (138->18), div. (22->4), fcn. (156->4), ass. (0->15)
	t42 = sin(pkin(10));
	t43 = cos(pkin(10));
	t44 = sin(qJ(1));
	t45 = cos(qJ(1));
	t37 = -t44 * t42 - t45 * t43;
	t35 = 0.1e1 / t37 ^ 2;
	t38 = t45 * t42 - t44 * t43;
	t50 = t38 ^ 2 * t35;
	t51 = t35 * t37;
	t32 = t38 * qJD(1);
	t34 = 0.1e1 / t37;
	t49 = t34 * t50;
	t48 = t32 * t51;
	t30 = 0.1e1 + t50;
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0.2e1 * (t34 * t37 + t50) / t30 ^ 2 * (t32 * t49 + t48) + (-0.2e1 * t48 + (t34 - 0.2e1 * t49 - t51) * t32) / t30, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:06:31
	% EndTime: 2019-10-10 00:06:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:06:32
	% EndTime: 2019-10-10 00:06:33
	% DurationCPUTime: 1.08s
	% Computational Cost: add. (1789->95), mult. (4580->216), div. (480->12), fcn. (5898->11), ass. (0->92)
	t218 = sin(pkin(10));
	t219 = cos(pkin(10));
	t220 = sin(qJ(1));
	t221 = cos(qJ(1));
	t149 = t221 * t218 - t220 * t219;
	t146 = t149 ^ 2;
	t162 = sin(qJ(4));
	t157 = t162 ^ 2;
	t164 = cos(qJ(4));
	t159 = 0.1e1 / t164 ^ 2;
	t199 = t157 * t159;
	t142 = t146 * t199 + 0.1e1;
	t148 = -t220 * t218 - t221 * t219;
	t144 = t148 * qJD(1);
	t156 = t162 * t157;
	t158 = 0.1e1 / t164;
	t198 = t158 * t162;
	t175 = qJD(4) * (t156 * t158 * t159 + t198);
	t211 = (t149 * t144 * t199 + t146 * t175) / t142 ^ 2;
	t223 = -0.2e1 * t211;
	t184 = 0.1e1 + t199;
	t222 = t149 * t184;
	t201 = t149 * t162;
	t139 = atan2(t201, t164);
	t137 = sin(t139);
	t138 = cos(t139);
	t128 = t137 * t201 + t138 * t164;
	t125 = 0.1e1 / t128;
	t161 = sin(qJ(5));
	t163 = cos(qJ(5));
	t196 = t163 * t164;
	t136 = -t148 * t196 + t149 * t161;
	t130 = 0.1e1 / t136;
	t126 = 0.1e1 / t128 ^ 2;
	t131 = 0.1e1 / t136 ^ 2;
	t147 = t148 ^ 2;
	t202 = t147 * t157;
	t119 = t126 * t202 + 0.1e1;
	t145 = t149 * qJD(1);
	t192 = qJD(4) * t164;
	t140 = 0.1e1 / t142;
	t186 = t149 * t192;
	t194 = qJD(4) * t149;
	t187 = t159 * t194;
	t203 = t144 * t162;
	t114 = ((t186 + t203) * t158 + t157 * t187) * t140;
	t204 = t138 * t162;
	t110 = (t114 * t149 - qJD(4)) * t204 + (t203 + (-t114 + t194) * t164) * t137;
	t216 = t110 * t125 * t126;
	t217 = (-t202 * t216 + (-t145 * t148 * t157 + t147 * t162 * t192) * t126) / t119 ^ 2;
	t193 = qJD(4) * t162;
	t172 = qJD(5) * t149 + t145 * t164 + t148 * t193;
	t191 = qJD(5) * t164;
	t178 = t148 * t191 + t144;
	t115 = t172 * t161 - t178 * t163;
	t197 = t161 * t164;
	t135 = -t148 * t197 - t149 * t163;
	t129 = t135 ^ 2;
	t124 = t129 * t131 + 0.1e1;
	t208 = t131 * t135;
	t116 = t178 * t161 + t172 * t163;
	t212 = t116 * t130 * t131;
	t215 = (t115 * t208 - t129 * t212) / t124 ^ 2;
	t123 = t140 * t222;
	t205 = t137 * t164;
	t112 = t149 * t205 - t204 + (t138 * t201 - t205) * t123;
	t214 = t112 * t126;
	t200 = t157 * t158;
	t188 = t149 * t200;
	t181 = t138 * t188;
	t206 = t137 * t162;
	t113 = (t206 + (t181 - t206) * t140) * t148;
	t213 = t113 * t126;
	t210 = t126 * t148;
	t209 = t130 * t161;
	t207 = t135 * t163;
	t195 = t123 - t149;
	t190 = -0.2e1 * t216;
	t189 = 0.2e1 * t215;
	t185 = t123 * t149 - 0.1e1;
	t183 = -0.2e1 * t162 * t217;
	t182 = 0.2e1 * t135 * t212;
	t177 = t149 * t191 + t145;
	t176 = t131 * t207 - t209;
	t174 = t176 * t162;
	t173 = qJD(5) * t148 + t144 * t164 - t149 * t193;
	t134 = t148 * t161 + t149 * t196;
	t133 = -t148 * t163 + t149 * t197;
	t121 = 0.1e1 / t124;
	t117 = 0.1e1 / t119;
	t109 = t222 * t223 + (t184 * t144 + 0.2e1 * t149 * t175) * t140;
	t1 = [t148 * t198 * t223 + (t184 * t148 * qJD(4) - t145 * t198) * t140, 0, 0, t109, 0, 0; t149 * t125 * t183 + (t125 * t186 + (t125 * t144 + (-t110 * t149 - t113 * t145) * t126) * t162) * t117 + (t183 * t213 + (t192 * t213 + (t113 * t190 + ((t137 * t192 + t181 * t223) * t148 + (-t137 * t145 + (t114 * t138 + 0.2e1 * t137 * t211) * t148) * t162 + (((-t114 * t188 - t192) * t148 + t145 * t162) * t137 + (-t145 * t188 + (t156 * t187 + t144 * t200 + (-t114 + 0.2e1 * t194) * t162) * t148) * t138) * t140) * t126) * t162) * t117) * t148, 0, 0, 0.2e1 * (t125 * t164 - t162 * t214) * t148 * t217 + ((t145 * t125 + (qJD(4) * t112 + t110) * t210) * t164 + (-t145 * t214 + (qJD(4) * t125 + t112 * t190 + ((t109 * t149 + t123 * t144) * t204 + (t195 * qJD(4) - t185 * t114) * t206) * t126) * t148 + ((-t109 + t144) * t137 + (t185 * qJD(4) - t195 * t114) * t138) * t164 * t210) * t162) * t117, 0, 0; (-t130 * t133 + t134 * t208) * t189 + (t134 * t182 + t177 * t130 * t163 + t173 * t209 + (t177 * t135 * t161 - t134 * t115 - t133 * t116 - t173 * t207) * t131) * t121, 0, 0, t148 * t174 * t189 + (t145 * t174 + (-t176 * t192 + ((qJD(5) * t130 + t182) * t163 + (-t115 * t163 + (qJD(5) * t135 - t116) * t161) * t131) * t162) * t148) * t121, -0.2e1 * t215 + 0.2e1 * (t115 * t131 * t121 + (-t121 * t212 - t131 * t215) * t135) * t135, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:06:32
	% EndTime: 2019-10-10 00:06:33
	% DurationCPUTime: 1.16s
	% Computational Cost: add. (2448->97), mult. (4991->216), div. (498->12), fcn. (6402->11), ass. (0->95)
	t249 = sin(pkin(10));
	t250 = cos(pkin(10));
	t251 = sin(qJ(1));
	t252 = cos(qJ(1));
	t178 = t252 * t249 - t251 * t250;
	t175 = t178 ^ 2;
	t194 = sin(qJ(4));
	t189 = t194 ^ 2;
	t195 = cos(qJ(4));
	t191 = 0.1e1 / t195 ^ 2;
	t227 = t189 * t191;
	t171 = t175 * t227 + 0.1e1;
	t177 = -t251 * t249 - t252 * t250;
	t173 = t177 * qJD(1);
	t188 = t194 * t189;
	t190 = 0.1e1 / t195;
	t226 = t190 * t194;
	t206 = qJD(4) * (t188 * t190 * t191 + t226);
	t242 = (t178 * t173 * t227 + t175 * t206) / t171 ^ 2;
	t254 = -0.2e1 * t242;
	t215 = 0.1e1 + t227;
	t253 = t178 * t215;
	t232 = t178 * t194;
	t168 = atan2(t232, t195);
	t166 = sin(t168);
	t167 = cos(t168);
	t157 = t166 * t232 + t167 * t195;
	t154 = 0.1e1 / t157;
	t193 = qJ(5) + qJ(6);
	t185 = sin(t193);
	t186 = cos(t193);
	t230 = t186 * t195;
	t165 = -t177 * t230 + t178 * t185;
	t159 = 0.1e1 / t165;
	t155 = 0.1e1 / t157 ^ 2;
	t160 = 0.1e1 / t165 ^ 2;
	t176 = t177 ^ 2;
	t233 = t176 * t189;
	t148 = t155 * t233 + 0.1e1;
	t174 = t178 * qJD(1);
	t222 = qJD(4) * t195;
	t169 = 0.1e1 / t171;
	t217 = t178 * t222;
	t224 = qJD(4) * t178;
	t218 = t191 * t224;
	t234 = t173 * t194;
	t145 = ((t217 + t234) * t190 + t189 * t218) * t169;
	t235 = t167 * t194;
	t140 = (t145 * t178 - qJD(4)) * t235 + (t234 + (-t145 + t224) * t195) * t166;
	t246 = t140 * t154 * t155;
	t248 = (-t233 * t246 + (-t174 * t177 * t189 + t176 * t194 * t222) * t155) / t148 ^ 2;
	t187 = qJD(5) + qJD(6);
	t223 = qJD(4) * t194;
	t203 = t174 * t195 + t177 * t223 + t178 * t187;
	t229 = t187 * t195;
	t209 = t177 * t229 + t173;
	t143 = t203 * t185 - t209 * t186;
	t231 = t185 * t195;
	t164 = -t177 * t231 - t178 * t186;
	t158 = t164 ^ 2;
	t151 = t158 * t160 + 0.1e1;
	t239 = t160 * t164;
	t144 = t209 * t185 + t203 * t186;
	t243 = t144 * t159 * t160;
	t247 = (t143 * t239 - t158 * t243) / t151 ^ 2;
	t153 = t169 * t253;
	t236 = t166 * t195;
	t141 = t178 * t236 - t235 + (t167 * t232 - t236) * t153;
	t245 = t141 * t155;
	t228 = t189 * t190;
	t219 = t178 * t228;
	t212 = t167 * t219;
	t237 = t166 * t194;
	t142 = (t237 + (t212 - t237) * t169) * t177;
	t244 = t142 * t155;
	t241 = t155 * t177;
	t240 = t159 * t185;
	t238 = t164 * t186;
	t225 = t153 - t178;
	t221 = 0.2e1 * t247;
	t220 = -0.2e1 * t246;
	t216 = t153 * t178 - 0.1e1;
	t214 = -0.2e1 * t194 * t248;
	t213 = 0.2e1 * t164 * t243;
	t208 = t178 * t229 + t174;
	t207 = t160 * t238 - t240;
	t205 = t207 * t194;
	t204 = t173 * t195 + t177 * t187 - t178 * t223;
	t163 = t177 * t185 + t178 * t230;
	t162 = -t177 * t186 + t178 * t231;
	t149 = 0.1e1 / t151;
	t146 = 0.1e1 / t148;
	t138 = t253 * t254 + (t215 * t173 + 0.2e1 * t178 * t206) * t169;
	t136 = -0.2e1 * t247 + 0.2e1 * (t143 * t160 * t149 + (-t149 * t243 - t160 * t247) * t164) * t164;
	t1 = [t177 * t226 * t254 + (t215 * t177 * qJD(4) - t174 * t226) * t169, 0, 0, t138, 0, 0; t178 * t154 * t214 + (t154 * t217 + (t154 * t173 + (-t140 * t178 - t142 * t174) * t155) * t194) * t146 + (t214 * t244 + (t222 * t244 + (t142 * t220 + ((t166 * t222 + t212 * t254) * t177 + (-t166 * t174 + (t145 * t167 + 0.2e1 * t166 * t242) * t177) * t194 + (((-t145 * t219 - t222) * t177 + t174 * t194) * t166 + (-t174 * t219 + (t188 * t218 + t173 * t228 + (-t145 + 0.2e1 * t224) * t194) * t177) * t167) * t169) * t155) * t194) * t146) * t177, 0, 0, 0.2e1 * (t154 * t195 - t194 * t245) * t177 * t248 + ((t174 * t154 + (qJD(4) * t141 + t140) * t241) * t195 + (-t174 * t245 + (qJD(4) * t154 + t141 * t220 + ((t138 * t178 + t153 * t173) * t235 + (t225 * qJD(4) - t216 * t145) * t237) * t155) * t177 + ((-t138 + t173) * t166 + (t216 * qJD(4) - t225 * t145) * t167) * t195 * t241) * t194) * t146, 0, 0; (-t159 * t162 + t163 * t239) * t221 + (t163 * t213 + t208 * t159 * t186 + t204 * t240 + (t208 * t164 * t185 - t163 * t143 - t162 * t144 - t204 * t238) * t160) * t149, 0, 0, t177 * t205 * t221 + (t174 * t205 + (-t207 * t222 + ((t159 * t187 + t213) * t186 + (-t143 * t186 + (t164 * t187 - t144) * t185) * t160) * t194) * t177) * t149, t136, t136;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end