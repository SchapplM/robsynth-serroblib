% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPPRRP4
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
%   Wie in S6RPPRRP4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:52
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRRP4_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP4_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRP4_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:52:56
	% EndTime: 2019-10-09 23:52:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:52:56
	% EndTime: 2019-10-09 23:52:56
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:52:56
	% EndTime: 2019-10-09 23:52:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:52:56
	% EndTime: 2019-10-09 23:52:56
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (52->9), mult. (138->18), div. (22->4), fcn. (156->4), ass. (0->15)
	t42 = sin(pkin(9));
	t43 = cos(pkin(9));
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
	% StartTime: 2019-10-09 23:52:56
	% EndTime: 2019-10-09 23:52:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:52:56
	% EndTime: 2019-10-09 23:52:57
	% DurationCPUTime: 1.11s
	% Computational Cost: add. (1789->95), mult. (4580->216), div. (480->12), fcn. (5898->11), ass. (0->92)
	t218 = sin(pkin(9));
	t219 = cos(pkin(9));
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
	% StartTime: 2019-10-09 23:52:56
	% EndTime: 2019-10-09 23:52:58
	% DurationCPUTime: 1.88s
	% Computational Cost: add. (4280->124), mult. (10807->281), div. (1114->15), fcn. (14382->11), ass. (0->112)
	t265 = sin(pkin(9));
	t266 = cos(pkin(9));
	t267 = sin(qJ(1));
	t268 = cos(qJ(1));
	t180 = -t267 * t265 - t268 * t266;
	t181 = t268 * t265 - t267 * t266;
	t196 = cos(qJ(5));
	t194 = sin(qJ(5));
	t197 = cos(qJ(4));
	t244 = t194 * t197;
	t164 = -t180 * t196 + t181 * t244;
	t195 = sin(qJ(4));
	t243 = t195 * t194;
	t155 = atan2(t164, -t243);
	t153 = sin(t155);
	t154 = cos(t155);
	t159 = t164 ^ 2;
	t188 = 0.1e1 / t194 ^ 2;
	t192 = 0.1e1 / t195 ^ 2;
	t246 = t188 * t192;
	t158 = t159 * t246 + 0.1e1;
	t156 = 0.1e1 / t158;
	t187 = 0.1e1 / t194;
	t225 = t187 * t192 * t197;
	t212 = -t164 * t225 - t181;
	t141 = t212 * t156;
	t241 = -t141 - t181;
	t272 = t241 * t195 * t153 - t154 * t197;
	t248 = t180 * t195;
	t242 = t196 * t197;
	t170 = -t180 * t242 + t181 * t194;
	t177 = t180 * qJD(1);
	t178 = t181 * qJD(1);
	t240 = qJD(4) * t195;
	t223 = t194 * t240;
	t148 = t170 * qJD(5) - t177 * t196 + t178 * t244 + t180 * t223;
	t169 = -t180 * t244 - t181 * t196;
	t191 = 0.1e1 / t195;
	t239 = qJD(4) * t197;
	t224 = t192 * t239;
	t237 = qJD(5) * t196;
	t247 = t187 * t191;
	t270 = -(t188 * t191 * t237 + t187 * t224) * t169 + t148 * t247;
	t257 = t153 * t164;
	t145 = -t154 * t243 + t257;
	t142 = 0.1e1 / t145;
	t161 = 0.1e1 / t170;
	t143 = 0.1e1 / t145 ^ 2;
	t162 = 0.1e1 / t170 ^ 2;
	t269 = -0.2e1 * t169;
	t160 = t169 ^ 2;
	t140 = t143 * t160 + 0.1e1;
	t259 = t148 * t143;
	t208 = -t194 * t239 - t195 * t237;
	t226 = t164 * t246;
	t165 = t180 * t194 + t181 * t242;
	t250 = t178 * t196;
	t146 = -t165 * qJD(5) - t177 * t244 + t181 * t223 - t250;
	t230 = t146 * t247;
	t135 = (-t208 * t226 + t230) * t156;
	t205 = -t135 * t164 - t208;
	t130 = (t135 * t243 - t146) * t153 - t205 * t154;
	t144 = t142 * t143;
	t263 = t130 * t144;
	t264 = (-t160 * t263 + t169 * t259) / t140 ^ 2;
	t179 = t180 ^ 2;
	t190 = t195 ^ 2;
	t249 = t179 * t190;
	t228 = t162 * t249;
	t152 = 0.1e1 + t228;
	t209 = -t178 * t197 - t180 * t240;
	t236 = qJD(5) * t197;
	t149 = (t180 * t236 + t177) * t194 + (qJD(5) * t181 - t209) * t196;
	t258 = t149 * t161 * t162;
	t216 = t249 * t258;
	t262 = (-t216 + (-t178 * t180 * t190 + t179 * t195 * t239) * t162) / t152 ^ 2;
	t189 = t187 * t188;
	t193 = t191 / t190;
	t221 = t192 * t237;
	t261 = (-t146 * t226 + (-t188 * t193 * t239 - t189 * t221) * t159) / t158 ^ 2;
	t260 = t143 * t169;
	t256 = t153 * t169;
	t254 = t154 * t164;
	t253 = t154 * t169;
	t251 = t162 * t180;
	t245 = t188 * t196;
	t238 = qJD(5) * t194;
	t235 = 0.2e1 * t264;
	t234 = 0.2e1 * t144 * t169;
	t233 = t195 * t262;
	t232 = t191 * t261;
	t231 = t143 * t256;
	t227 = t164 * t247;
	t222 = t196 * t239;
	t220 = -0.2e1 * t142 * t264;
	t219 = t143 * t235;
	t218 = t187 * t232;
	t217 = t248 * t258;
	t213 = t181 * t236 + t178;
	t211 = -t164 * t245 + t165 * t187;
	t210 = -t178 * t195 + t180 * t239;
	t207 = -qJD(5) * t180 - t177 * t197 + t181 * t240;
	t147 = t213 * t194 + t207 * t196;
	t150 = 0.1e1 / t152;
	t138 = 0.1e1 / t140;
	t137 = t211 * t191 * t156;
	t133 = (-t153 + (t154 * t227 + t153) * t156) * t169;
	t132 = -t141 * t254 + t272 * t194;
	t131 = -t154 * t195 * t196 + t153 * t165 - (t153 * t243 + t254) * t137;
	t129 = 0.2e1 * t212 * t261 + (-t146 * t225 + t177 - (t188 * t197 * t221 + (0.2e1 * t193 * t197 ^ 2 + t191) * t187 * qJD(4)) * t164) * t156;
	t127 = 0.2e1 * t211 * t232 + (t211 * t224 + (-t146 * t245 + t147 * t187 + (t165 * t245 - (0.2e1 * t189 * t196 ^ 2 + t187) * t164) * qJD(5)) * t191) * t156;
	t1 = [t270 * t156 + t218 * t269, 0, 0, t129, t127, 0; t164 * t220 + ((-t207 * t194 + t213 * t196) * t142 + (-t130 * t164 - t133 * t148) * t143) * t138 + (t133 * t219 + (0.2e1 * t133 * t263 + (-t148 * t156 + t148 - (-t135 * t156 * t227 - 0.2e1 * t261) * t169) * t143 * t153 + (-(-0.2e1 * t164 * t218 - t135) * t260 + (-(t135 - t230) * t169 - t270 * t164) * t143 * t156) * t154) * t138) * t169, 0, 0, t132 * t169 * t219 + (-(t129 * t254 - (-t135 * t257 - t146 * t154) * t141) * t260 + (t130 * t234 - t259) * t132 + (t142 * t248 - t272 * t260) * t237) * t138 + (t220 * t248 + ((t180 * qJD(4) * t142 - (t241 * qJD(4) + t135) * t231) * t197 + (-t178 * t142 + (-t180 * t130 - (t129 - t177) * t256 - (t241 * t135 + qJD(4)) * t253) * t143) * t195) * t138) * t194, (t131 * t260 - t142 * t170) * t235 + (-t131 * t259 + t149 * t142 + (t131 * t234 - t143 * t170) * t130 - (-t222 + t195 * t238 + t127 * t164 + t137 * t146 + (-t137 * t243 + t165) * t135) * t143 * t253 - (-t147 + (t127 * t194 + t135 * t196) * t195 - t205 * t137) * t231) * t138, 0; 0.2e1 * (t161 * t181 + t165 * t251) * t233 + (0.2e1 * t165 * t217 + (-t177 * t195 - t181 * t239) * t161 + ((t147 * t180 + t181 * t149) * t195 - t210 * t165) * t162) * t150, 0, 0, 0.2e1 * (-t161 * t180 * t197 + t196 * t228) * t262 + (0.2e1 * t196 * t216 + t209 * t161 + ((-t149 * t197 + 0.2e1 * t190 * t250) * t180 + (t190 * t238 - 0.2e1 * t195 * t222) * t179) * t162) * t150, t233 * t251 * t269 + (t217 * t269 + (t148 * t248 + t210 * t169) * t162) * t150, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end