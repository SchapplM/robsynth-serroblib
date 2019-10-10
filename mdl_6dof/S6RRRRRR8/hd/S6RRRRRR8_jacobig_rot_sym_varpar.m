% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRRR8
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:29
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRRRR8_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR8_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_jacobig_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:29:30
	% EndTime: 2019-10-10 13:29:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:29:30
	% EndTime: 2019-10-10 13:29:30
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:29:30
	% EndTime: 2019-10-10 13:29:30
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t43 = sin(pkin(6));
	t1 = [0, sin(qJ(1)) * t43, 0, 0, 0, 0; 0, -cos(qJ(1)) * t43, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:29:31
	% EndTime: 2019-10-10 13:29:31
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (6->6), mult. (17->14), div. (0->0), fcn. (28->8), ass. (0->12)
	t94 = cos(pkin(6));
	t97 = cos(qJ(2));
	t101 = t94 * t97;
	t92 = sin(pkin(6));
	t96 = sin(qJ(1));
	t100 = t96 * t92;
	t98 = cos(qJ(1));
	t99 = t98 * t92;
	t95 = sin(qJ(2));
	t93 = cos(pkin(7));
	t91 = sin(pkin(7));
	t1 = [0, t100, -(-t96 * t101 - t98 * t95) * t91 + t93 * t100, 0, 0, 0; 0, -t99, -(t98 * t101 - t96 * t95) * t91 - t93 * t99, 0, 0, 0; 1, t94, -t92 * t97 * t91 + t94 * t93, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:29:31
	% EndTime: 2019-10-10 13:29:31
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (16->14), mult. (48->31), div. (0->0), fcn. (72->10), ass. (0->19)
	t140 = sin(pkin(6));
	t145 = sin(qJ(1));
	t154 = t145 * t140;
	t144 = sin(qJ(2));
	t153 = t145 * t144;
	t147 = cos(qJ(2));
	t152 = t145 * t147;
	t148 = cos(qJ(1));
	t151 = t148 * t140;
	t150 = t148 * t144;
	t149 = t148 * t147;
	t146 = cos(qJ(3));
	t143 = sin(qJ(3));
	t142 = cos(pkin(6));
	t141 = cos(pkin(7));
	t139 = sin(pkin(7));
	t138 = -t142 * t152 - t150;
	t137 = t142 * t149 - t153;
	t1 = [0, t154, -t138 * t139 + t141 * t154, (-t142 * t153 + t149) * t143 + (-t138 * t141 - t139 * t154) * t146, 0, 0; 0, -t151, -t137 * t139 - t141 * t151, (t142 * t150 + t152) * t143 + (-t137 * t141 + t139 * t151) * t146, 0, 0; 1, t142, -t140 * t147 * t139 + t142 * t141, -t142 * t139 * t146 + (-t141 * t146 * t147 + t143 * t144) * t140, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:29:31
	% EndTime: 2019-10-10 13:29:31
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (26->14), mult. (79->31), div. (0->0), fcn. (116->10), ass. (0->22)
	t158 = sin(pkin(6));
	t163 = sin(qJ(1));
	t172 = t163 * t158;
	t162 = sin(qJ(2));
	t171 = t163 * t162;
	t165 = cos(qJ(2));
	t170 = t163 * t165;
	t166 = cos(qJ(1));
	t169 = t166 * t158;
	t168 = t166 * t162;
	t167 = t166 * t165;
	t164 = cos(qJ(3));
	t161 = sin(qJ(3));
	t160 = cos(pkin(6));
	t159 = cos(pkin(7));
	t157 = sin(pkin(7));
	t156 = -t160 * t170 - t168;
	t155 = t160 * t167 - t171;
	t154 = -t160 * t157 * t164 + (-t159 * t164 * t165 + t161 * t162) * t158;
	t153 = (-t160 * t171 + t167) * t161 + (-t156 * t159 - t157 * t172) * t164;
	t152 = (t160 * t168 + t170) * t161 + (-t155 * t159 + t157 * t169) * t164;
	t1 = [0, t172, -t156 * t157 + t159 * t172, t153, t153, 0; 0, -t169, -t155 * t157 - t159 * t169, t152, t152, 0; 1, t160, -t158 * t165 * t157 + t160 * t159, t154, t154, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:29:32
	% EndTime: 2019-10-10 13:29:32
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (50->22), mult. (131->45), div. (0->0), fcn. (189->12), ass. (0->34)
	t220 = sin(pkin(7));
	t223 = cos(pkin(6));
	t239 = t220 * t223;
	t222 = cos(pkin(7));
	t228 = cos(qJ(2));
	t238 = t222 * t228;
	t221 = sin(pkin(6));
	t226 = sin(qJ(1));
	t237 = t226 * t221;
	t225 = sin(qJ(2));
	t236 = t226 * t225;
	t235 = t226 * t228;
	t229 = cos(qJ(1));
	t234 = t229 * t221;
	t233 = t229 * t225;
	t232 = t229 * t228;
	t213 = t223 * t232 - t236;
	t231 = -t213 * t222 + t220 * t234;
	t215 = -t223 * t235 - t233;
	t230 = t215 * t222 + t220 * t237;
	t227 = cos(qJ(3));
	t224 = sin(qJ(3));
	t219 = qJ(4) + qJ(5);
	t218 = cos(t219);
	t217 = sin(t219);
	t216 = -t223 * t236 + t232;
	t214 = t223 * t233 + t235;
	t212 = -t221 * t228 * t220 + t223 * t222;
	t211 = -t215 * t220 + t222 * t237;
	t210 = -t213 * t220 - t222 * t234;
	t209 = -t227 * t239 + (t224 * t225 - t227 * t238) * t221;
	t208 = t216 * t224 - t230 * t227;
	t207 = t214 * t224 + t231 * t227;
	t1 = [0, t237, t211, t208, t208, (t216 * t227 + t230 * t224) * t217 - t211 * t218; 0, -t234, t210, t207, t207, (t214 * t227 - t231 * t224) * t217 - t210 * t218; 1, t223, t212, t209, t209, (t224 * t239 + (t224 * t238 + t225 * t227) * t221) * t217 - t212 * t218;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end