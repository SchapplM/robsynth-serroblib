% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S4PRRR7
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S4PRRR7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% JaD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 12:35
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S4PRRR7_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_jacobiaD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR7_jacobiaD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PRRR7_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:34:58
	% EndTime: 2019-12-29 12:34:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:34:58
	% EndTime: 2019-12-29 12:34:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:34:58
	% EndTime: 2019-12-29 12:34:58
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (46->7), mult. (159->21), div. (18->4), fcn. (175->5), ass. (0->15)
	t39 = cos(pkin(8));
	t41 = sin(qJ(2));
	t42 = cos(qJ(2));
	t45 = sin(pkin(8)) * cos(pkin(4));
	t37 = t39 * t42 - t41 * t45;
	t34 = 0.1e1 / t37 ^ 2;
	t49 = qJD(2) * t34;
	t36 = t39 * t41 + t42 * t45;
	t33 = t36 ^ 2;
	t30 = t33 * t34 + 0.1e1;
	t46 = t37 * t49;
	t47 = t36 / t37 * t49;
	t48 = (t33 * t47 + t36 * t46) / t30 ^ 2;
	t28 = 0.1e1 / t30;
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, -0.2e1 * t48 + 0.2e1 * (t28 * t46 + (t28 * t47 - t34 * t48) * t36) * t36, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:34:59
	% EndTime: 2019-12-29 12:34:59
	% DurationCPUTime: 0.83s
	% Computational Cost: add. (756->55), mult. (2271->133), div. (423->14), fcn. (2956->11), ass. (0->65)
	t138 = sin(qJ(2));
	t140 = cos(qJ(2));
	t135 = sin(pkin(8));
	t165 = cos(pkin(4));
	t154 = t135 * t165;
	t164 = cos(pkin(8));
	t126 = -t138 * t154 + t164 * t140;
	t137 = sin(qJ(3));
	t139 = cos(qJ(3));
	t136 = sin(pkin(4));
	t159 = t135 * t136;
	t148 = -t126 * t137 + t139 * t159;
	t169 = t148 * qJD(3);
	t151 = t165 * t164;
	t122 = t135 * t138 - t140 * t151;
	t158 = t136 * t140;
	t112 = atan2(-t122, -t158);
	t110 = sin(t112);
	t111 = cos(t112);
	t97 = -t110 * t122 - t111 * t158;
	t94 = 0.1e1 / t97;
	t109 = t126 * t139 + t137 * t159;
	t105 = 0.1e1 / t109;
	t132 = 0.1e1 / t140;
	t106 = 0.1e1 / t109 ^ 2;
	t133 = 0.1e1 / t140 ^ 2;
	t95 = 0.1e1 / t97 ^ 2;
	t104 = t148 ^ 2;
	t101 = t104 * t106 + 0.1e1;
	t147 = -t164 * t138 - t140 * t154;
	t118 = t147 * qJD(2);
	t102 = t109 * qJD(3) + t118 * t137;
	t162 = t106 * t148;
	t103 = t118 * t139 + t169;
	t163 = t103 * t105 * t106;
	t168 = 0.1e1 / t101 ^ 2 * (-t102 * t162 - t104 * t163);
	t124 = t135 * t140 + t138 * t151;
	t160 = t133 * t138;
	t155 = t122 * t160;
	t149 = t124 * t132 + t155;
	t120 = t122 ^ 2;
	t131 = 0.1e1 / t136 ^ 2;
	t115 = t120 * t131 * t133 + 0.1e1;
	t113 = 0.1e1 / t115;
	t130 = 0.1e1 / t136;
	t161 = t113 * t130;
	t90 = t149 * t161;
	t167 = t122 * t90;
	t166 = t147 * t95;
	t157 = qJD(2) * t138;
	t156 = -0.2e1 * t168;
	t150 = -t105 * t137 - t139 * t162;
	t134 = t132 * t133;
	t121 = t147 ^ 2;
	t119 = t126 * qJD(2);
	t117 = t124 * qJD(2);
	t116 = t122 * qJD(2);
	t99 = 0.1e1 / t101;
	t96 = t94 * t95;
	t93 = t121 * t95 + 0.1e1;
	t89 = (qJD(2) * t155 + t117 * t132) * t161;
	t87 = (t136 * t138 - t167) * t111 + (t90 * t158 - t124) * t110;
	t86 = (-t122 * t89 + t136 * t157) * t111 + (t89 * t158 - t117) * t110;
	t85 = (-0.2e1 * t149 * (t117 * t122 * t133 + t120 * t134 * t157) * t131 / t115 ^ 2 + (t117 * t160 - t116 * t132 + (t124 * t160 + (0.2e1 * t134 * t138 ^ 2 + t132) * t122) * qJD(2)) * t113) * t130;
	t1 = [0, t85, 0, 0; 0, 0.2e1 * (-t126 * t94 - t87 * t166) / t93 ^ 2 * (-t121 * t96 * t86 - t119 * t166) + (-t87 * t119 * t95 + t118 * t94 + (-0.2e1 * t147 * t87 * t96 - t126 * t95) * t86 - (-(-t117 * t90 - t122 * t85 - t124 * t89 + (t89 * t90 + qJD(2)) * t158) * t111 - (t89 * t167 + t116 + (t140 * t85 + (-qJD(2) * t90 - t89) * t138) * t136) * t110) * t166) / t93, 0, 0; 0, t150 * t99 * t119 - (t150 * t156 + ((-qJD(3) * t105 + 0.2e1 * t148 * t163) * t139 + (t102 * t139 + (t103 + t169) * t137) * t106) * t99) * t147, t156 - 0.2e1 * (t102 * t106 * t99 - (-t106 * t168 - t99 * t163) * t148) * t148, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:34:59
	% EndTime: 2019-12-29 12:35:01
	% DurationCPUTime: 2.07s
	% Computational Cost: add. (3002->109), mult. (9085->226), div. (559->12), fcn. (11668->13), ass. (0->102)
	t207 = sin(pkin(8));
	t209 = cos(pkin(8));
	t213 = sin(qJ(2));
	t210 = cos(pkin(4));
	t216 = cos(qJ(2));
	t242 = t210 * t216;
	t197 = -t207 * t213 + t209 * t242;
	t190 = t197 * qJD(2);
	t243 = t210 * t213;
	t198 = t207 * t216 + t209 * t243;
	t212 = sin(qJ(3));
	t208 = sin(pkin(4));
	t246 = t208 * t212;
	t232 = t209 * t246;
	t215 = cos(qJ(3));
	t239 = qJD(3) * t215;
	t162 = -qJD(3) * t232 + t190 * t212 + t198 * t239;
	t245 = t208 * t215;
	t182 = t198 * t212 + t209 * t245;
	t180 = t182 ^ 2;
	t201 = -t210 * t215 + t213 * t246;
	t195 = 0.1e1 / t201 ^ 2;
	t176 = t180 * t195 + 0.1e1;
	t174 = 0.1e1 / t176;
	t202 = t210 * t212 + t213 * t245;
	t240 = qJD(2) * t216;
	t231 = t208 * t240;
	t187 = t202 * qJD(3) + t212 * t231;
	t194 = 0.1e1 / t201;
	t250 = t182 * t195;
	t146 = (-t162 * t194 + t187 * t250) * t174;
	t177 = atan2(-t182, t201);
	t172 = sin(t177);
	t173 = cos(t177);
	t228 = -t172 * t201 - t173 * t182;
	t142 = t228 * t146 - t172 * t162 + t173 * t187;
	t156 = -t172 * t182 + t173 * t201;
	t153 = 0.1e1 / t156;
	t154 = 0.1e1 / t156 ^ 2;
	t262 = t142 * t153 * t154;
	t233 = t207 * t243;
	t200 = t209 * t216 - t233;
	t225 = -t200 * t212 + t207 * t245;
	t261 = -0.2e1 * t225 * t262;
	t244 = t208 * t216;
	t224 = -t194 * t197 + t244 * t250;
	t260 = t212 * t224;
	t249 = t187 * t194 * t195;
	t259 = -0.2e1 * (t162 * t250 - t180 * t249) / t176 ^ 2;
	t186 = t200 * t215 + t207 * t246;
	t199 = t207 * t242 + t209 * t213;
	t211 = sin(qJ(4));
	t214 = cos(qJ(4));
	t171 = t186 * t214 + t199 * t211;
	t167 = 0.1e1 / t171;
	t168 = 0.1e1 / t171 ^ 2;
	t192 = t199 * qJD(2);
	t165 = t225 * qJD(3) - t192 * t215;
	t193 = -qJD(2) * t233 + t209 * t240;
	t157 = t171 * qJD(4) + t165 * t211 - t193 * t214;
	t170 = t186 * t211 - t199 * t214;
	t166 = t170 ^ 2;
	t161 = t166 * t168 + 0.1e1;
	t254 = t168 * t170;
	t238 = qJD(4) * t170;
	t158 = t165 * t214 + t193 * t211 - t238;
	t256 = t158 * t167 * t168;
	t258 = (t157 * t254 - t166 * t256) / t161 ^ 2;
	t257 = t154 * t225;
	t255 = t167 * t211;
	t253 = t170 * t214;
	t252 = t172 * t225;
	t251 = t173 * t225;
	t248 = t199 * t212;
	t247 = t199 * t215;
	t241 = qJD(2) * t213;
	t181 = t225 ^ 2;
	t152 = t181 * t154 + 0.1e1;
	t164 = t186 * qJD(3) - t192 * t212;
	t237 = 0.2e1 * (-t164 * t257 - t181 * t262) / t152 ^ 2;
	t236 = -0.2e1 * t258;
	t234 = t170 * t256;
	t230 = -0.2e1 * t182 * t249;
	t229 = qJD(4) * t247 - t192;
	t227 = t168 * t253 - t255;
	t184 = t198 * t215 - t232;
	t226 = -t184 * t194 + t202 * t250;
	t223 = qJD(3) * t248 + qJD(4) * t200 - t193 * t215;
	t191 = t198 * qJD(2);
	t188 = -t201 * qJD(3) + t215 * t231;
	t179 = t200 * t211 - t214 * t247;
	t178 = -t200 * t214 - t211 * t247;
	t163 = -t182 * qJD(3) + t190 * t215;
	t159 = 0.1e1 / t161;
	t149 = 0.1e1 / t152;
	t148 = t174 * t260;
	t147 = t226 * t174;
	t144 = (-t172 * t197 + t173 * t244) * t212 + t228 * t148;
	t143 = t228 * t147 - t172 * t184 + t173 * t202;
	t141 = t226 * t259 + (t202 * t230 - t163 * t194 + (t162 * t202 + t182 * t188 + t184 * t187) * t195) * t174;
	t139 = t259 * t260 + (t224 * t239 + (t230 * t244 + t191 * t194 + (t187 * t197 + (t162 * t216 - t182 * t241) * t208) * t195) * t212) * t174;
	t1 = [0, t139, t141, 0; 0, (-t144 * t257 + t153 * t248) * t237 + ((-t193 * t212 - t199 * t239) * t153 + t144 * t261 + (-t144 * t164 + t248 * t142 + (-t139 * t182 - t148 * t162 + (-t212 * t241 + t216 * t239) * t208 + (-t148 * t201 - t197 * t212) * t146) * t251 + (-t197 * t239 - t139 * t201 - t148 * t187 + t191 * t212 + (t148 * t182 - t212 * t244) * t146) * t252) * t154) * t149, (-t143 * t257 - t153 * t186) * t237 + (t143 * t261 + t165 * t153 + (-t186 * t142 - t143 * t164 + (-t141 * t182 - t147 * t162 + t188 + (-t147 * t201 - t184) * t146) * t251 + (-t141 * t201 - t147 * t187 - t163 + (t147 * t182 - t202) * t146) * t252) * t154) * t149, 0; 0, 0.2e1 * (-t167 * t178 + t179 * t254) * t258 + (0.2e1 * t179 * t234 - t229 * t167 * t214 + t223 * t255 + (-t229 * t170 * t211 - t179 * t157 - t178 * t158 - t223 * t253) * t168) * t159, -t227 * t225 * t236 + (t227 * t164 - ((-qJD(4) * t167 - 0.2e1 * t234) * t214 + (t157 * t214 + (t158 - t238) * t211) * t168) * t225) * t159, t236 + 0.2e1 * (t157 * t168 * t159 + (-t159 * t256 - t168 * t258) * t170) * t170;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,4);
end