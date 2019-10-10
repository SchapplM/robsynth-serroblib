% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:41
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRPR11_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR11_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_jacobiR_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:41:20
	% EndTime: 2019-10-10 01:41:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:41:20
	% EndTime: 2019-10-10 01:41:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0, 0; t9, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0, 0; -t8, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:41:20
	% EndTime: 2019-10-10 01:41:20
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (4->4), mult. (14->10), div. (0->0), fcn. (24->6), ass. (0->11)
	t38 = sin(pkin(12));
	t42 = sin(qJ(1));
	t47 = t42 * t38;
	t40 = cos(pkin(12));
	t46 = t42 * t40;
	t43 = cos(qJ(1));
	t45 = t43 * t38;
	t44 = t43 * t40;
	t41 = cos(pkin(6));
	t39 = sin(pkin(6));
	t1 = [-t41 * t45 - t46, 0, 0, 0, 0, 0; -t41 * t47 + t44, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t41 * t44 + t47, 0, 0, 0, 0, 0; -t41 * t46 - t45, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t43 * t39, 0, 0, 0, 0, 0; t42 * t39, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:41:20
	% EndTime: 2019-10-10 01:41:20
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (40->17), mult. (122->36), div. (0->0), fcn. (174->10), ass. (0->27)
	t103 = cos(qJ(1));
	t96 = sin(pkin(6));
	t111 = t103 * t96;
	t101 = sin(qJ(1));
	t99 = cos(pkin(6));
	t110 = t103 * t99;
	t94 = sin(pkin(12));
	t97 = cos(pkin(12));
	t88 = t101 * t94 - t97 * t110;
	t95 = sin(pkin(7));
	t98 = cos(pkin(7));
	t117 = t95 * t111 + t88 * t98;
	t100 = sin(qJ(3));
	t102 = cos(qJ(3));
	t89 = t101 * t97 + t94 * t110;
	t116 = t89 * t100 + t117 * t102;
	t114 = t94 * t96;
	t113 = t101 * t96;
	t112 = t101 * t99;
	t107 = t96 * t97 * t98 + t95 * t99;
	t90 = -t103 * t94 - t97 * t112;
	t106 = t95 * t113 + t90 * t98;
	t104 = t117 * t100 - t89 * t102;
	t91 = t103 * t97 - t94 * t112;
	t87 = t106 * t100 + t91 * t102;
	t86 = -t91 * t100 + t106 * t102;
	t1 = [t104, 0, t86, 0, 0, 0; t87, 0, -t116, 0, 0, 0; 0, 0, -t100 * t114 + t107 * t102, 0, 0, 0; t116, 0, -t87, 0, 0, 0; t86, 0, t104, 0, 0, 0; 0, 0, -t107 * t100 - t102 * t114, 0, 0, 0; t98 * t111 - t88 * t95, 0, 0, 0, 0, 0; t98 * t113 - t90 * t95, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:41:20
	% EndTime: 2019-10-10 01:41:21
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (105->30), mult. (307->59), div. (0->0), fcn. (430->12), ass. (0->40)
	t160 = cos(pkin(6));
	t155 = sin(pkin(12));
	t166 = cos(qJ(1));
	t170 = t166 * t155;
	t158 = cos(pkin(12));
	t163 = sin(qJ(1));
	t171 = t163 * t158;
	t150 = t160 * t170 + t171;
	t162 = sin(qJ(3));
	t165 = cos(qJ(3));
	t169 = t166 * t158;
	t172 = t163 * t155;
	t149 = -t160 * t169 + t172;
	t156 = sin(pkin(7));
	t159 = cos(pkin(7));
	t157 = sin(pkin(6));
	t174 = t157 * t166;
	t167 = t149 * t159 + t156 * t174;
	t139 = -t150 * t165 + t167 * t162;
	t144 = -t149 * t156 + t159 * t174;
	t161 = sin(qJ(4));
	t164 = cos(qJ(4));
	t181 = t139 * t161 - t144 * t164;
	t180 = t139 * t164 + t144 * t161;
	t176 = t156 * t160;
	t175 = t157 * t163;
	t173 = t159 * t165;
	t168 = t156 * t175;
	t137 = -t150 * t162 - t167 * t165;
	t152 = -t160 * t172 + t169;
	t151 = -t160 * t171 - t170;
	t148 = -t157 * t158 * t156 + t160 * t159;
	t146 = -t151 * t156 + t159 * t175;
	t143 = t162 * t176 + (t158 * t159 * t162 + t155 * t165) * t157;
	t142 = t165 * t176 + (-t155 * t162 + t158 * t173) * t157;
	t141 = t152 * t165 + (t151 * t159 + t168) * t162;
	t140 = -t151 * t173 + t152 * t162 - t165 * t168;
	t136 = t141 * t164 + t146 * t161;
	t135 = -t141 * t161 + t146 * t164;
	t1 = [t180, 0, -t140 * t164, t135, 0, 0; t136, 0, t137 * t164, t181, 0, 0; 0, 0, t142 * t164, -t143 * t161 + t148 * t164, 0, 0; -t181, 0, t140 * t161, -t136, 0, 0; t135, 0, -t137 * t161, t180, 0, 0; 0, 0, -t142 * t161, -t143 * t164 - t148 * t161, 0, 0; t137, 0, t141, 0, 0, 0; t140, 0, -t139, 0, 0, 0; 0, 0, t143, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:41:21
	% EndTime: 2019-10-10 01:41:21
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (175->38), mult. (514->82), div. (0->0), fcn. (710->14), ass. (0->45)
	t195 = cos(pkin(6));
	t189 = sin(pkin(12));
	t201 = cos(qJ(1));
	t209 = t201 * t189;
	t193 = cos(pkin(12));
	t198 = sin(qJ(1));
	t210 = t198 * t193;
	t184 = t195 * t209 + t210;
	t197 = sin(qJ(3));
	t200 = cos(qJ(3));
	t208 = t201 * t193;
	t211 = t198 * t189;
	t183 = -t195 * t208 + t211;
	t190 = sin(pkin(7));
	t194 = cos(pkin(7));
	t191 = sin(pkin(6));
	t214 = t191 * t201;
	t206 = t183 * t194 + t190 * t214;
	t173 = -t184 * t200 + t206 * t197;
	t179 = -t183 * t190 + t194 * t214;
	t196 = sin(qJ(4));
	t199 = cos(qJ(4));
	t165 = t173 * t196 - t179 * t199;
	t166 = t173 * t199 + t179 * t196;
	t205 = t195 * t210 + t209;
	t215 = t191 * t198;
	t219 = -t190 * t215 + t205 * t194;
	t188 = sin(pkin(13));
	t217 = t188 * t199;
	t216 = t190 * t195;
	t192 = cos(pkin(13));
	t213 = t192 * t199;
	t212 = t193 * t194;
	t203 = -t184 * t197 - t206 * t200;
	t202 = t205 * t190 + t194 * t215;
	t185 = -t195 * t211 + t208;
	t182 = -t191 * t193 * t190 + t195 * t194;
	t178 = t197 * t216 + (t189 * t200 + t197 * t212) * t191;
	t177 = t200 * t216 + (-t189 * t197 + t200 * t212) * t191;
	t175 = t185 * t200 - t219 * t197;
	t174 = t185 * t197 + t219 * t200;
	t169 = -t178 * t196 + t182 * t199;
	t168 = t175 * t199 + t202 * t196;
	t167 = t175 * t196 - t202 * t199;
	t1 = [t166 * t192 + t188 * t203, 0, -t174 * t213 + t175 * t188, -t167 * t192, 0, 0; t168 * t192 + t174 * t188, 0, -t173 * t188 + t203 * t213, t165 * t192, 0, 0; 0, 0, t177 * t213 + t178 * t188, t169 * t192, 0, 0; -t166 * t188 + t192 * t203, 0, t174 * t217 + t175 * t192, t167 * t188, 0, 0; -t168 * t188 + t174 * t192, 0, -t173 * t192 - t203 * t217, -t165 * t188, 0, 0; 0, 0, -t177 * t217 + t178 * t192, -t169 * t188, 0, 0; t165, 0, -t174 * t196, t168, 0, 0; t167, 0, t203 * t196, -t166, 0, 0; 0, 0, t177 * t196, t178 * t199 + t182 * t196, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:41:21
	% EndTime: 2019-10-10 01:41:21
	% DurationCPUTime: 0.35s
	% Computational Cost: add. (275->48), mult. (692->91), div. (0->0), fcn. (956->14), ass. (0->55)
	t233 = cos(pkin(6));
	t228 = sin(pkin(12));
	t239 = cos(qJ(1));
	t247 = t239 * t228;
	t231 = cos(pkin(12));
	t236 = sin(qJ(1));
	t248 = t236 * t231;
	t220 = t233 * t247 + t248;
	t235 = sin(qJ(3));
	t238 = cos(qJ(3));
	t246 = t239 * t231;
	t249 = t236 * t228;
	t219 = -t233 * t246 + t249;
	t232 = cos(pkin(7));
	t229 = sin(pkin(7));
	t230 = sin(pkin(6));
	t251 = t230 * t239;
	t245 = t229 * t251;
	t243 = t219 * t232 + t245;
	t208 = -t220 * t238 + t243 * t235;
	t214 = -t219 * t229 + t232 * t251;
	t234 = sin(qJ(4));
	t237 = cos(qJ(4));
	t200 = t208 * t237 + t214 * t234;
	t227 = pkin(13) + qJ(6);
	t225 = sin(t227);
	t263 = t200 * t225;
	t226 = cos(t227);
	t262 = t200 * t226;
	t198 = t208 * t234 - t214 * t237;
	t242 = t233 * t248 + t247;
	t252 = t230 * t236;
	t259 = -t229 * t252 + t242 * t232;
	t258 = t220 * t235;
	t256 = t225 * t237;
	t255 = t226 * t237;
	t254 = t229 * t233;
	t253 = t230 * t231;
	t250 = t232 * t238;
	t240 = t242 * t229 + t232 * t252;
	t221 = -t233 * t249 + t246;
	t218 = -t229 * t253 + t233 * t232;
	t213 = t235 * t254 + (t231 * t232 * t235 + t228 * t238) * t230;
	t212 = t230 * t228 * t235 - t238 * t254 - t250 * t253;
	t210 = t221 * t238 - t259 * t235;
	t209 = t221 * t235 + t259 * t238;
	t207 = -t243 * t238 - t258;
	t205 = t219 * t250 + t238 * t245 + t258;
	t204 = t213 * t237 + t218 * t234;
	t203 = -t213 * t234 + t218 * t237;
	t202 = t210 * t237 + t240 * t234;
	t201 = t210 * t234 - t240 * t237;
	t197 = t202 * t226 + t209 * t225;
	t196 = -t202 * t225 + t209 * t226;
	t1 = [t207 * t225 + t262, 0, -t209 * t255 + t210 * t225, -t201 * t226, 0, t196; t197, 0, -t205 * t255 - t208 * t225, t198 * t226, 0, t205 * t226 + t263; 0, 0, -t212 * t255 + t213 * t225, t203 * t226, 0, -t204 * t225 + t212 * t226; t207 * t226 - t263, 0, t209 * t256 + t210 * t226, t201 * t225, 0, -t197; t196, 0, t205 * t256 - t208 * t226, -t198 * t225, 0, -t205 * t225 + t262; 0, 0, t212 * t256 + t213 * t226, -t203 * t225, 0, -t204 * t226 - t212 * t225; t198, 0, -t209 * t234, t202, 0, 0; t201, 0, -t205 * t234, -t200, 0, 0; 0, 0, -t212 * t234, t204, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end