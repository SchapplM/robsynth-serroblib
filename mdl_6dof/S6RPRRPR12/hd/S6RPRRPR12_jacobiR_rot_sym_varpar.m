% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRRPR12
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:43
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRPR12_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR12_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_jacobiR_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:43:25
	% EndTime: 2019-10-10 01:43:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:43:25
	% EndTime: 2019-10-10 01:43:25
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
	% StartTime: 2019-10-10 01:43:25
	% EndTime: 2019-10-10 01:43:25
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
	% StartTime: 2019-10-10 01:43:25
	% EndTime: 2019-10-10 01:43:26
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
	% StartTime: 2019-10-10 01:43:26
	% EndTime: 2019-10-10 01:43:26
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
	% StartTime: 2019-10-10 01:43:26
	% EndTime: 2019-10-10 01:43:26
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (105->30), mult. (307->59), div. (0->0), fcn. (430->12), ass. (0->40)
	t180 = cos(pkin(6));
	t175 = sin(pkin(12));
	t186 = cos(qJ(1));
	t190 = t186 * t175;
	t178 = cos(pkin(12));
	t183 = sin(qJ(1));
	t191 = t183 * t178;
	t170 = t180 * t190 + t191;
	t182 = sin(qJ(3));
	t185 = cos(qJ(3));
	t189 = t186 * t178;
	t192 = t183 * t175;
	t169 = -t180 * t189 + t192;
	t176 = sin(pkin(7));
	t179 = cos(pkin(7));
	t177 = sin(pkin(6));
	t194 = t177 * t186;
	t187 = t169 * t179 + t176 * t194;
	t159 = -t170 * t185 + t187 * t182;
	t164 = -t169 * t176 + t179 * t194;
	t181 = sin(qJ(4));
	t184 = cos(qJ(4));
	t201 = t159 * t181 - t164 * t184;
	t200 = -t159 * t184 - t164 * t181;
	t196 = t176 * t180;
	t195 = t177 * t183;
	t193 = t179 * t185;
	t188 = t176 * t195;
	t157 = -t170 * t182 - t187 * t185;
	t172 = -t180 * t192 + t189;
	t171 = -t180 * t191 - t190;
	t168 = -t177 * t178 * t176 + t180 * t179;
	t166 = -t171 * t176 + t179 * t195;
	t163 = t182 * t196 + (t178 * t179 * t182 + t175 * t185) * t177;
	t162 = t185 * t196 + (-t175 * t182 + t178 * t193) * t177;
	t161 = t172 * t185 + (t171 * t179 + t188) * t182;
	t160 = -t171 * t193 + t172 * t182 - t185 * t188;
	t156 = t161 * t184 + t166 * t181;
	t155 = t161 * t181 - t166 * t184;
	t1 = [t157, 0, t161, 0, 0, 0; t160, 0, -t159, 0, 0, 0; 0, 0, t163, 0, 0, 0; t200, 0, t160 * t184, t155, 0, 0; -t156, 0, -t157 * t184, -t201, 0, 0; 0, 0, -t162 * t184, t163 * t181 - t168 * t184, 0, 0; t201, 0, -t160 * t181, t156, 0, 0; t155, 0, t157 * t181, t200, 0, 0; 0, 0, t162 * t181, t163 * t184 + t168 * t181, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:43:26
	% EndTime: 2019-10-10 01:43:26
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (234->47), mult. (692->91), div. (0->0), fcn. (956->14), ass. (0->54)
	t226 = cos(pkin(6));
	t221 = sin(pkin(12));
	t234 = cos(qJ(1));
	t241 = t234 * t221;
	t224 = cos(pkin(12));
	t230 = sin(qJ(1));
	t242 = t230 * t224;
	t216 = t226 * t241 + t242;
	t229 = sin(qJ(3));
	t233 = cos(qJ(3));
	t240 = t234 * t224;
	t243 = t230 * t221;
	t215 = -t226 * t240 + t243;
	t225 = cos(pkin(7));
	t222 = sin(pkin(7));
	t223 = sin(pkin(6));
	t247 = t223 * t234;
	t239 = t222 * t247;
	t237 = t215 * t225 + t239;
	t204 = -t216 * t233 + t237 * t229;
	t209 = -t215 * t222 + t225 * t247;
	t228 = sin(qJ(4));
	t232 = cos(qJ(4));
	t196 = t204 * t228 - t209 * t232;
	t227 = sin(qJ(6));
	t258 = t196 * t227;
	t231 = cos(qJ(6));
	t257 = t196 * t231;
	t256 = t204 * t232 + t209 * t228;
	t236 = t226 * t242 + t241;
	t248 = t223 * t230;
	t253 = -t222 * t248 + t236 * t225;
	t252 = t216 * t229;
	t250 = t222 * t226;
	t249 = t223 * t224;
	t246 = t225 * t233;
	t245 = t227 * t228;
	t244 = t228 * t231;
	t217 = -t226 * t243 + t240;
	t214 = -t222 * t249 + t226 * t225;
	t211 = t236 * t222 + t225 * t248;
	t208 = t229 * t250 + (t224 * t225 * t229 + t221 * t233) * t223;
	t207 = t223 * t221 * t229 - t233 * t250 - t246 * t249;
	t206 = t217 * t233 - t253 * t229;
	t205 = t217 * t229 + t253 * t233;
	t203 = -t237 * t233 - t252;
	t201 = t215 * t246 + t233 * t239 + t252;
	t200 = t208 * t232 + t214 * t228;
	t199 = t208 * t228 - t214 * t232;
	t198 = t206 * t232 + t211 * t228;
	t197 = t206 * t228 - t211 * t232;
	t193 = t197 * t227 + t205 * t231;
	t192 = t197 * t231 - t205 * t227;
	t1 = [t203 * t231 + t258, 0, -t205 * t245 + t206 * t231, t198 * t227, 0, t192; t193, 0, -t201 * t245 - t204 * t231, -t256 * t227, 0, -t201 * t227 - t257; 0, 0, -t207 * t245 + t208 * t231, t200 * t227, 0, t199 * t231 - t207 * t227; -t203 * t227 + t257, 0, -t205 * t244 - t206 * t227, t198 * t231, 0, -t193; t192, 0, -t201 * t244 + t204 * t227, -t256 * t231, 0, -t201 * t231 + t258; 0, 0, -t207 * t244 - t208 * t227, t200 * t231, 0, -t199 * t227 - t207 * t231; t256, 0, -t205 * t232, -t197, 0, 0; t198, 0, -t201 * t232, t196, 0, 0; 0, 0, -t207 * t232, -t199, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end