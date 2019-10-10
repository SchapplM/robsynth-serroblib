% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:07
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPRR13_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR13_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_jacobiR_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:07:44
	% EndTime: 2019-10-10 01:07:44
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:07:44
	% EndTime: 2019-10-10 01:07:44
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
	% StartTime: 2019-10-10 01:07:44
	% EndTime: 2019-10-10 01:07:44
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
	% StartTime: 2019-10-10 01:07:44
	% EndTime: 2019-10-10 01:07:44
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
	% StartTime: 2019-10-10 01:07:44
	% EndTime: 2019-10-10 01:07:44
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (40->18), mult. (122->36), div. (0->0), fcn. (174->10), ass. (0->29)
	t124 = cos(pkin(6));
	t119 = sin(pkin(12));
	t128 = cos(qJ(1));
	t133 = t128 * t119;
	t122 = cos(pkin(12));
	t126 = sin(qJ(1));
	t134 = t126 * t122;
	t115 = t124 * t133 + t134;
	t125 = sin(qJ(3));
	t127 = cos(qJ(3));
	t132 = t128 * t122;
	t135 = t126 * t119;
	t114 = -t124 * t132 + t135;
	t120 = sin(pkin(7));
	t123 = cos(pkin(7));
	t121 = sin(pkin(6));
	t137 = t121 * t128;
	t131 = t114 * t123 + t120 * t137;
	t141 = t115 * t125 + t131 * t127;
	t139 = t120 * t124;
	t138 = t121 * t126;
	t136 = t122 * t123;
	t116 = -t124 * t134 - t133;
	t130 = t116 * t123 + t120 * t138;
	t129 = t115 * t127 - t131 * t125;
	t117 = -t124 * t135 + t132;
	t113 = t117 * t127 + t130 * t125;
	t112 = t117 * t125 - t130 * t127;
	t1 = [-t114 * t120 + t123 * t137, 0, 0, 0, 0, 0; -t116 * t120 + t123 * t138, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t129, 0, t112, 0, 0, 0; -t113, 0, t141, 0, 0, 0; 0, 0, -t127 * t139 + (t119 * t125 - t127 * t136) * t121, 0, 0, 0; -t141, 0, t113, 0, 0, 0; t112, 0, t129, 0, 0, 0; 0, 0, t125 * t139 + (t119 * t127 + t125 * t136) * t121, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:07:45
	% EndTime: 2019-10-10 01:07:45
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (102->30), mult. (307->57), div. (0->0), fcn. (430->12), ass. (0->40)
	t158 = cos(pkin(6));
	t153 = sin(pkin(12));
	t164 = cos(qJ(1));
	t168 = t164 * t153;
	t156 = cos(pkin(12));
	t161 = sin(qJ(1));
	t169 = t161 * t156;
	t149 = t158 * t168 + t169;
	t160 = sin(qJ(3));
	t163 = cos(qJ(3));
	t167 = t164 * t156;
	t170 = t161 * t153;
	t148 = -t158 * t167 + t170;
	t154 = sin(pkin(7));
	t157 = cos(pkin(7));
	t155 = sin(pkin(6));
	t172 = t155 * t164;
	t166 = t148 * t157 + t154 * t172;
	t137 = t149 * t160 + t166 * t163;
	t144 = -t148 * t154 + t157 * t172;
	t159 = sin(qJ(5));
	t162 = cos(qJ(5));
	t181 = -t137 * t159 + t144 * t162;
	t180 = t137 * t162 + t144 * t159;
	t177 = -t149 * t163 + t166 * t160;
	t174 = t154 * t158;
	t173 = t155 * t161;
	t171 = t156 * t157;
	t150 = -t158 * t169 - t168;
	t165 = t150 * t157 + t154 * t173;
	t151 = -t158 * t170 + t167;
	t147 = -t155 * t156 * t154 + t158 * t157;
	t146 = -t150 * t154 + t157 * t173;
	t143 = t160 * t174 + (t153 * t163 + t160 * t171) * t155;
	t142 = -t163 * t174 + (t153 * t160 - t163 * t171) * t155;
	t141 = t151 * t163 + t165 * t160;
	t140 = t151 * t160 - t165 * t163;
	t136 = t140 * t159 + t146 * t162;
	t135 = t140 * t162 - t146 * t159;
	t1 = [t181, 0, t141 * t159, 0, t135, 0; t136, 0, -t177 * t159, 0, t180, 0; 0, 0, t143 * t159, 0, t142 * t162 - t147 * t159, 0; -t180, 0, t141 * t162, 0, -t136, 0; t135, 0, -t177 * t162, 0, t181, 0; 0, 0, t143 * t162, 0, -t142 * t159 - t147 * t162, 0; t177, 0, -t140, 0, 0, 0; t141, 0, -t137, 0, 0, 0; 0, 0, -t142, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:07:45
	% EndTime: 2019-10-10 01:07:45
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (240->48), mult. (692->94), div. (0->0), fcn. (956->14), ass. (0->57)
	t219 = cos(pkin(6));
	t214 = sin(pkin(12));
	t227 = cos(qJ(1));
	t232 = t227 * t214;
	t217 = cos(pkin(12));
	t223 = sin(qJ(1));
	t233 = t223 * t217;
	t207 = t219 * t232 + t233;
	t222 = sin(qJ(3));
	t226 = cos(qJ(3));
	t231 = t227 * t217;
	t234 = t223 * t214;
	t206 = -t219 * t231 + t234;
	t218 = cos(pkin(7));
	t215 = sin(pkin(7));
	t216 = sin(pkin(6));
	t238 = t216 * t227;
	t229 = t215 * t238;
	t228 = t206 * t218 + t229;
	t195 = -t207 * t226 + t228 * t222;
	t220 = sin(qJ(6));
	t247 = t195 * t220;
	t224 = cos(qJ(6));
	t246 = t195 * t224;
	t200 = -t206 * t215 + t218 * t238;
	t221 = sin(qJ(5));
	t245 = t200 * t221;
	t225 = cos(qJ(5));
	t244 = t200 * t225;
	t243 = t207 * t222;
	t241 = t215 * t219;
	t240 = t216 * t217;
	t239 = t216 * t223;
	t237 = t218 * t226;
	t236 = t220 * t221;
	t235 = t221 * t224;
	t230 = t215 * t239;
	t209 = -t219 * t234 + t231;
	t208 = -t219 * t233 - t232;
	t205 = -t215 * t240 + t219 * t218;
	t202 = -t208 * t215 + t218 * t239;
	t199 = t222 * t241 + (t217 * t218 * t222 + t214 * t226) * t216;
	t198 = t216 * t214 * t222 - t226 * t241 - t237 * t240;
	t197 = t209 * t226 + (t208 * t218 + t230) * t222;
	t196 = -t208 * t237 + t209 * t222 - t226 * t230;
	t194 = -t228 * t226 - t243;
	t192 = t206 * t237 + t226 * t229 + t243;
	t191 = t198 * t221 + t205 * t225;
	t190 = t198 * t225 - t205 * t221;
	t188 = t196 * t221 + t202 * t225;
	t187 = -t196 * t225 + t202 * t221;
	t186 = t192 * t221 - t244;
	t185 = t192 * t225 + t245;
	t184 = t194 * t221 + t244;
	t183 = t188 * t224 + t197 * t220;
	t182 = -t188 * t220 + t197 * t224;
	t1 = [t184 * t224 + t247, 0, -t196 * t220 + t197 * t235, 0, -t187 * t224, t182; t183, 0, -t192 * t220 - t195 * t235, 0, t185 * t224, -t186 * t220 - t246; 0, 0, -t198 * t220 + t199 * t235, 0, t190 * t224, -t191 * t220 + t199 * t224; -t184 * t220 + t246, 0, -t196 * t224 - t197 * t236, 0, t187 * t220, -t183; t182, 0, -t192 * t224 + t195 * t236, 0, -t185 * t220, -t186 * t224 + t247; 0, 0, -t198 * t224 - t199 * t236, 0, -t190 * t220, -t191 * t224 - t199 * t220; -t194 * t225 + t245, 0, -t197 * t225, 0, t188, 0; t187, 0, t195 * t225, 0, t186, 0; 0, 0, -t199 * t225, 0, t191, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end