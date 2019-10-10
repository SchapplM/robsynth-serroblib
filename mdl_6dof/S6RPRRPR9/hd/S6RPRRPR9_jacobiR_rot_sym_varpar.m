% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRRPR9
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
% Datum: 2019-10-10 01:37
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRPR9_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR9_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_jacobiR_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:37:31
	% EndTime: 2019-10-10 01:37:31
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:37:31
	% EndTime: 2019-10-10 01:37:31
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
	% StartTime: 2019-10-10 01:37:31
	% EndTime: 2019-10-10 01:37:31
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
	% StartTime: 2019-10-10 01:37:32
	% EndTime: 2019-10-10 01:37:32
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
	% StartTime: 2019-10-10 01:37:32
	% EndTime: 2019-10-10 01:37:32
	% DurationCPUTime: 0.17s
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
	% StartTime: 2019-10-10 01:37:32
	% EndTime: 2019-10-10 01:37:32
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (131->31), mult. (307->59), div. (0->0), fcn. (430->12), ass. (0->41)
	t173 = cos(pkin(6));
	t168 = sin(pkin(12));
	t177 = cos(qJ(1));
	t181 = t177 * t168;
	t171 = cos(pkin(12));
	t175 = sin(qJ(1));
	t182 = t175 * t171;
	t160 = t173 * t181 + t182;
	t174 = sin(qJ(3));
	t176 = cos(qJ(3));
	t180 = t177 * t171;
	t183 = t175 * t168;
	t159 = -t173 * t180 + t183;
	t169 = sin(pkin(7));
	t172 = cos(pkin(7));
	t170 = sin(pkin(6));
	t185 = t170 * t177;
	t178 = t159 * t172 + t169 * t185;
	t149 = -t160 * t176 + t178 * t174;
	t154 = -t159 * t169 + t172 * t185;
	t167 = qJ(4) + pkin(13);
	t165 = sin(t167);
	t166 = cos(t167);
	t192 = t149 * t165 - t154 * t166;
	t191 = t149 * t166 + t154 * t165;
	t187 = t169 * t173;
	t186 = t170 * t175;
	t184 = t172 * t176;
	t179 = t169 * t186;
	t147 = -t160 * t174 - t178 * t176;
	t162 = -t173 * t183 + t180;
	t161 = -t173 * t182 - t181;
	t158 = -t170 * t171 * t169 + t173 * t172;
	t156 = -t161 * t169 + t172 * t186;
	t153 = t174 * t187 + (t171 * t172 * t174 + t168 * t176) * t170;
	t152 = t176 * t187 + (-t168 * t174 + t171 * t184) * t170;
	t151 = t162 * t176 + (t161 * t172 + t179) * t174;
	t150 = -t161 * t184 + t162 * t174 - t176 * t179;
	t146 = t151 * t166 + t156 * t165;
	t145 = -t151 * t165 + t156 * t166;
	t1 = [t191, 0, -t150 * t166, t145, 0, 0; t146, 0, t147 * t166, t192, 0, 0; 0, 0, t152 * t166, -t153 * t165 + t158 * t166, 0, 0; -t192, 0, t150 * t165, -t146, 0, 0; t145, 0, -t147 * t165, t191, 0, 0; 0, 0, -t152 * t165, -t153 * t166 - t158 * t165, 0, 0; t147, 0, t151, 0, 0, 0; t150, 0, -t149, 0, 0, 0; 0, 0, t153, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:37:32
	% EndTime: 2019-10-10 01:37:33
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (288->48), mult. (692->91), div. (0->0), fcn. (956->14), ass. (0->55)
	t229 = cos(pkin(6));
	t224 = sin(pkin(12));
	t235 = cos(qJ(1));
	t243 = t235 * t224;
	t227 = cos(pkin(12));
	t232 = sin(qJ(1));
	t244 = t232 * t227;
	t216 = t229 * t243 + t244;
	t231 = sin(qJ(3));
	t234 = cos(qJ(3));
	t242 = t235 * t227;
	t245 = t232 * t224;
	t215 = -t229 * t242 + t245;
	t228 = cos(pkin(7));
	t225 = sin(pkin(7));
	t226 = sin(pkin(6));
	t247 = t226 * t235;
	t241 = t225 * t247;
	t239 = t215 * t228 + t241;
	t204 = -t216 * t234 + t239 * t231;
	t210 = -t215 * t225 + t228 * t247;
	t223 = qJ(4) + pkin(13);
	t221 = sin(t223);
	t222 = cos(t223);
	t196 = t204 * t222 + t210 * t221;
	t230 = sin(qJ(6));
	t259 = t196 * t230;
	t233 = cos(qJ(6));
	t258 = t196 * t233;
	t194 = t204 * t221 - t210 * t222;
	t238 = t229 * t244 + t243;
	t248 = t226 * t232;
	t255 = -t225 * t248 + t238 * t228;
	t254 = t216 * t231;
	t252 = t222 * t230;
	t251 = t222 * t233;
	t250 = t225 * t229;
	t249 = t226 * t227;
	t246 = t228 * t234;
	t236 = t238 * t225 + t228 * t248;
	t217 = -t229 * t245 + t242;
	t214 = -t225 * t249 + t229 * t228;
	t209 = t231 * t250 + (t227 * t228 * t231 + t224 * t234) * t226;
	t208 = t226 * t224 * t231 - t234 * t250 - t246 * t249;
	t206 = t217 * t234 - t255 * t231;
	t205 = t217 * t231 + t255 * t234;
	t203 = -t239 * t234 - t254;
	t201 = t215 * t246 + t234 * t241 + t254;
	t200 = t209 * t222 + t214 * t221;
	t199 = -t209 * t221 + t214 * t222;
	t198 = t206 * t222 + t236 * t221;
	t197 = t206 * t221 - t236 * t222;
	t193 = t198 * t233 + t205 * t230;
	t192 = -t198 * t230 + t205 * t233;
	t1 = [t203 * t230 + t258, 0, -t205 * t251 + t206 * t230, -t197 * t233, 0, t192; t193, 0, -t201 * t251 - t204 * t230, t194 * t233, 0, t201 * t233 + t259; 0, 0, -t208 * t251 + t209 * t230, t199 * t233, 0, -t200 * t230 + t208 * t233; t203 * t233 - t259, 0, t205 * t252 + t206 * t233, t197 * t230, 0, -t193; t192, 0, t201 * t252 - t204 * t233, -t194 * t230, 0, -t201 * t230 + t258; 0, 0, t208 * t252 + t209 * t233, -t199 * t230, 0, -t200 * t233 - t208 * t230; t194, 0, -t205 * t221, t198, 0, 0; t197, 0, -t201 * t221, -t196, 0, 0; 0, 0, -t208 * t221, t200, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end