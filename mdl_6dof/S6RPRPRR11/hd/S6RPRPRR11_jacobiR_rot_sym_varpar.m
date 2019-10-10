% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:04
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPRR11_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR11_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_jacobiR_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:04:00
	% EndTime: 2019-10-10 01:04:00
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:04:00
	% EndTime: 2019-10-10 01:04:00
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
	% StartTime: 2019-10-10 01:04:00
	% EndTime: 2019-10-10 01:04:00
	% DurationCPUTime: 0.04s
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
	% StartTime: 2019-10-10 01:04:00
	% EndTime: 2019-10-10 01:04:00
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
	% StartTime: 2019-10-10 01:04:00
	% EndTime: 2019-10-10 01:04:00
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (69->25), mult. (203->52), div. (0->0), fcn. (284->12), ass. (0->34)
	t143 = cos(pkin(6));
	t137 = sin(pkin(12));
	t147 = cos(qJ(1));
	t151 = t147 * t137;
	t141 = cos(pkin(12));
	t145 = sin(qJ(1));
	t152 = t145 * t141;
	t131 = t143 * t151 + t152;
	t144 = sin(qJ(3));
	t146 = cos(qJ(3));
	t150 = t147 * t141;
	t153 = t145 * t137;
	t130 = -t143 * t150 + t153;
	t138 = sin(pkin(7));
	t142 = cos(pkin(7));
	t139 = sin(pkin(6));
	t155 = t139 * t147;
	t148 = t130 * t142 + t138 * t155;
	t123 = -t131 * t146 + t148 * t144;
	t157 = t138 * t143;
	t156 = t139 * t145;
	t154 = t142 * t146;
	t149 = t138 * t156;
	t122 = -t131 * t144 - t148 * t146;
	t140 = cos(pkin(13));
	t136 = sin(pkin(13));
	t133 = -t143 * t153 + t150;
	t132 = -t143 * t152 - t151;
	t128 = -t132 * t138 + t142 * t156;
	t127 = -t130 * t138 + t142 * t155;
	t126 = t146 * t157 + (-t137 * t144 + t141 * t154) * t139;
	t125 = t133 * t146 + (t132 * t142 + t149) * t144;
	t124 = -t132 * t154 + t133 * t144 - t146 * t149;
	t1 = [t123 * t140 + t127 * t136, 0, -t124 * t140, 0, 0, 0; t125 * t140 + t128 * t136, 0, t122 * t140, 0, 0, 0; 0, 0, t126 * t140, 0, 0, 0; -t123 * t136 + t127 * t140, 0, t124 * t136, 0, 0, 0; -t125 * t136 + t128 * t140, 0, -t122 * t136, 0, 0, 0; 0, 0, -t126 * t136, 0, 0, 0; t122, 0, t125, 0, 0, 0; t124, 0, -t123, 0, 0, 0; 0, 0, t144 * t157 + (t141 * t142 * t144 + t137 * t146) * t139, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:04:00
	% EndTime: 2019-10-10 01:04:00
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (131->31), mult. (307->59), div. (0->0), fcn. (430->12), ass. (0->41)
	t170 = cos(pkin(6));
	t165 = sin(pkin(12));
	t174 = cos(qJ(1));
	t178 = t174 * t165;
	t168 = cos(pkin(12));
	t172 = sin(qJ(1));
	t179 = t172 * t168;
	t157 = t170 * t178 + t179;
	t171 = sin(qJ(3));
	t173 = cos(qJ(3));
	t177 = t174 * t168;
	t180 = t172 * t165;
	t156 = -t170 * t177 + t180;
	t166 = sin(pkin(7));
	t169 = cos(pkin(7));
	t167 = sin(pkin(6));
	t182 = t167 * t174;
	t175 = t156 * t169 + t166 * t182;
	t146 = -t157 * t173 + t175 * t171;
	t151 = -t156 * t166 + t169 * t182;
	t164 = pkin(13) + qJ(5);
	t162 = sin(t164);
	t163 = cos(t164);
	t189 = t146 * t162 - t151 * t163;
	t188 = t146 * t163 + t151 * t162;
	t184 = t166 * t170;
	t183 = t167 * t172;
	t181 = t169 * t173;
	t176 = t166 * t183;
	t144 = -t157 * t171 - t175 * t173;
	t159 = -t170 * t180 + t177;
	t158 = -t170 * t179 - t178;
	t155 = -t167 * t168 * t166 + t170 * t169;
	t153 = -t158 * t166 + t169 * t183;
	t150 = t171 * t184 + (t168 * t169 * t171 + t165 * t173) * t167;
	t149 = t173 * t184 + (-t165 * t171 + t168 * t181) * t167;
	t148 = t159 * t173 + (t158 * t169 + t176) * t171;
	t147 = -t158 * t181 + t159 * t171 - t173 * t176;
	t143 = t148 * t163 + t153 * t162;
	t142 = -t148 * t162 + t153 * t163;
	t1 = [t188, 0, -t147 * t163, 0, t142, 0; t143, 0, t144 * t163, 0, t189, 0; 0, 0, t149 * t163, 0, -t150 * t162 + t155 * t163, 0; -t189, 0, t147 * t162, 0, -t143, 0; t142, 0, -t144 * t162, 0, t188, 0; 0, 0, -t149 * t162, 0, -t150 * t163 - t155 * t162, 0; t144, 0, t148, 0, 0, 0; t147, 0, -t146, 0, 0, 0; 0, 0, t150, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:04:01
	% EndTime: 2019-10-10 01:04:01
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (288->48), mult. (692->91), div. (0->0), fcn. (956->14), ass. (0->55)
	t230 = cos(pkin(6));
	t225 = sin(pkin(12));
	t236 = cos(qJ(1));
	t244 = t236 * t225;
	t228 = cos(pkin(12));
	t233 = sin(qJ(1));
	t245 = t233 * t228;
	t217 = t230 * t244 + t245;
	t232 = sin(qJ(3));
	t235 = cos(qJ(3));
	t243 = t236 * t228;
	t246 = t233 * t225;
	t216 = -t230 * t243 + t246;
	t229 = cos(pkin(7));
	t226 = sin(pkin(7));
	t227 = sin(pkin(6));
	t248 = t227 * t236;
	t242 = t226 * t248;
	t240 = t216 * t229 + t242;
	t205 = -t217 * t235 + t240 * t232;
	t211 = -t216 * t226 + t229 * t248;
	t224 = pkin(13) + qJ(5);
	t222 = sin(t224);
	t223 = cos(t224);
	t197 = t205 * t223 + t211 * t222;
	t231 = sin(qJ(6));
	t260 = t197 * t231;
	t234 = cos(qJ(6));
	t259 = t197 * t234;
	t195 = t205 * t222 - t211 * t223;
	t239 = t230 * t245 + t244;
	t249 = t227 * t233;
	t256 = -t226 * t249 + t239 * t229;
	t255 = t217 * t232;
	t253 = t223 * t231;
	t252 = t223 * t234;
	t251 = t226 * t230;
	t250 = t227 * t228;
	t247 = t229 * t235;
	t237 = t239 * t226 + t229 * t249;
	t218 = -t230 * t246 + t243;
	t215 = -t226 * t250 + t230 * t229;
	t210 = t232 * t251 + (t228 * t229 * t232 + t225 * t235) * t227;
	t209 = t227 * t225 * t232 - t235 * t251 - t247 * t250;
	t207 = t218 * t235 - t256 * t232;
	t206 = t218 * t232 + t256 * t235;
	t204 = -t240 * t235 - t255;
	t202 = t216 * t247 + t235 * t242 + t255;
	t201 = t210 * t223 + t215 * t222;
	t200 = -t210 * t222 + t215 * t223;
	t199 = t207 * t223 + t237 * t222;
	t198 = t207 * t222 - t237 * t223;
	t194 = t199 * t234 + t206 * t231;
	t193 = -t199 * t231 + t206 * t234;
	t1 = [t204 * t231 + t259, 0, -t206 * t252 + t207 * t231, 0, -t198 * t234, t193; t194, 0, -t202 * t252 - t205 * t231, 0, t195 * t234, t202 * t234 + t260; 0, 0, -t209 * t252 + t210 * t231, 0, t200 * t234, -t201 * t231 + t209 * t234; t204 * t234 - t260, 0, t206 * t253 + t207 * t234, 0, t198 * t231, -t194; t193, 0, t202 * t253 - t205 * t234, 0, -t195 * t231, -t202 * t231 + t259; 0, 0, t209 * t253 + t210 * t234, 0, -t200 * t231, -t201 * t234 - t209 * t231; t195, 0, -t206 * t222, 0, t199, 0; t198, 0, -t202 * t222, 0, -t197, 0; 0, 0, -t209 * t222, 0, t201, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end