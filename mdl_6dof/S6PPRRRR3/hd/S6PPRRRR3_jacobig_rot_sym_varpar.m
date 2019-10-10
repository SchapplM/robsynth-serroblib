% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PPRRRR3
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:22
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PPRRRR3_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRRR3_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_jacobig_rot_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:20
	% EndTime: 2019-10-09 21:22:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:20
	% EndTime: 2019-10-09 21:22:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:20
	% EndTime: 2019-10-09 21:22:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:21
	% EndTime: 2019-10-09 21:22:21
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->5), mult. (15->13), div. (0->0), fcn. (23->8), ass. (0->11)
	t41 = sin(pkin(6));
	t44 = cos(pkin(7));
	t47 = t41 * t44;
	t42 = cos(pkin(14));
	t45 = cos(pkin(6));
	t46 = t42 * t45;
	t43 = cos(pkin(13));
	t40 = sin(pkin(7));
	t39 = sin(pkin(13));
	t38 = sin(pkin(14));
	t1 = [0, 0, -(-t43 * t38 - t39 * t46) * t40 + t39 * t47, 0, 0, 0; 0, 0, -(-t39 * t38 + t43 * t46) * t40 - t43 * t47, 0, 0, 0; 0, 0, -t41 * t42 * t40 + t45 * t44, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:21
	% EndTime: 2019-10-09 21:22:21
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (23->16), mult. (67->38), div. (0->0), fcn. (96->12), ass. (0->22)
	t102 = sin(pkin(13));
	t110 = cos(pkin(6));
	t116 = t102 * t110;
	t104 = sin(pkin(7));
	t105 = sin(pkin(6));
	t115 = t104 * t105;
	t109 = cos(pkin(7));
	t114 = t105 * t109;
	t107 = cos(pkin(13));
	t113 = t107 * t110;
	t112 = cos(qJ(3));
	t111 = sin(qJ(3));
	t108 = cos(pkin(8));
	t106 = cos(pkin(14));
	t103 = sin(pkin(8));
	t101 = sin(pkin(14));
	t100 = -t107 * t101 - t106 * t116;
	t99 = -t102 * t101 + t106 * t113;
	t98 = -t106 * t115 + t110 * t109;
	t97 = -t100 * t104 + t102 * t114;
	t96 = -t99 * t104 - t107 * t114;
	t1 = [0, 0, t97, -(-(-t101 * t116 + t107 * t106) * t111 + (t100 * t109 + t102 * t115) * t112) * t103 + t97 * t108, 0, 0; 0, 0, t96, -(-(t101 * t113 + t102 * t106) * t111 + (-t107 * t115 + t99 * t109) * t112) * t103 + t96 * t108, 0, 0; 0, 0, t98, -(t110 * t104 * t112 + (t106 * t109 * t112 - t101 * t111) * t105) * t103 + t98 * t108, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:21
	% EndTime: 2019-10-09 21:22:21
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (54->26), mult. (159->58), div. (0->0), fcn. (222->14), ass. (0->33)
	t157 = sin(pkin(13));
	t165 = cos(pkin(6));
	t177 = t157 * t165;
	t159 = sin(pkin(7));
	t160 = sin(pkin(6));
	t176 = t159 * t160;
	t175 = t159 * t165;
	t164 = cos(pkin(7));
	t174 = t160 * t164;
	t161 = cos(pkin(14));
	t173 = t161 * t164;
	t162 = cos(pkin(13));
	t172 = t162 * t165;
	t156 = sin(pkin(14));
	t152 = -t157 * t156 + t161 * t172;
	t171 = t152 * t164 - t162 * t176;
	t154 = -t162 * t156 - t161 * t177;
	t170 = t154 * t164 + t157 * t176;
	t169 = cos(qJ(3));
	t168 = cos(qJ(4));
	t167 = sin(qJ(3));
	t166 = sin(qJ(4));
	t163 = cos(pkin(8));
	t158 = sin(pkin(8));
	t155 = -t156 * t177 + t162 * t161;
	t153 = t156 * t172 + t157 * t161;
	t151 = -t161 * t176 + t165 * t164;
	t150 = -t154 * t159 + t157 * t174;
	t149 = -t152 * t159 - t162 * t174;
	t148 = t169 * t175 + (-t156 * t167 + t169 * t173) * t160;
	t147 = -t155 * t167 + t170 * t169;
	t146 = -t153 * t167 + t171 * t169;
	t1 = [0, 0, t150, -t147 * t158 + t150 * t163, (t155 * t169 + t170 * t167) * t166 + (-t147 * t163 - t150 * t158) * t168, 0; 0, 0, t149, -t146 * t158 + t149 * t163, (t153 * t169 + t171 * t167) * t166 + (-t146 * t163 - t149 * t158) * t168, 0; 0, 0, t151, -t148 * t158 + t151 * t163, (t167 * t175 + (t156 * t169 + t167 * t173) * t160) * t166 + (-t148 * t163 - t151 * t158) * t168, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:22
	% EndTime: 2019-10-09 21:22:22
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (106->32), mult. (309->70), div. (0->0), fcn. (427->16), ass. (0->44)
	t206 = sin(pkin(13));
	t214 = cos(pkin(6));
	t231 = t206 * t214;
	t208 = sin(pkin(7));
	t209 = sin(pkin(6));
	t230 = t208 * t209;
	t229 = t208 * t214;
	t213 = cos(pkin(7));
	t228 = t209 * t213;
	t210 = cos(pkin(14));
	t227 = t210 * t213;
	t211 = cos(pkin(13));
	t226 = t211 * t214;
	t205 = sin(pkin(14));
	t202 = t205 * t226 + t206 * t210;
	t217 = sin(qJ(3));
	t220 = cos(qJ(3));
	t201 = -t206 * t205 + t210 * t226;
	t222 = t201 * t213 - t211 * t230;
	t191 = -t202 * t217 + t222 * t220;
	t198 = -t201 * t208 - t211 * t228;
	t207 = sin(pkin(8));
	t212 = cos(pkin(8));
	t225 = t191 * t212 + t198 * t207;
	t204 = -t205 * t231 + t211 * t210;
	t203 = -t211 * t205 - t210 * t231;
	t221 = t203 * t213 + t206 * t230;
	t193 = -t204 * t217 + t221 * t220;
	t199 = -t203 * t208 + t206 * t228;
	t224 = t193 * t212 + t199 * t207;
	t196 = t220 * t229 + (-t205 * t217 + t220 * t227) * t209;
	t200 = -t210 * t230 + t214 * t213;
	t223 = t196 * t212 + t200 * t207;
	t219 = cos(qJ(4));
	t218 = cos(qJ(5));
	t216 = sin(qJ(4));
	t215 = sin(qJ(5));
	t197 = t217 * t229 + (t205 * t220 + t217 * t227) * t209;
	t195 = -t196 * t207 + t200 * t212;
	t194 = t204 * t220 + t221 * t217;
	t192 = t202 * t220 + t222 * t217;
	t190 = -t193 * t207 + t199 * t212;
	t189 = -t191 * t207 + t198 * t212;
	t1 = [0, 0, t199, t190, t194 * t216 - t224 * t219, (t194 * t219 + t224 * t216) * t215 - t190 * t218; 0, 0, t198, t189, t192 * t216 - t225 * t219, (t192 * t219 + t225 * t216) * t215 - t189 * t218; 0, 0, t200, t195, t197 * t216 - t223 * t219, (t197 * t219 + t223 * t216) * t215 - t195 * t218;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end