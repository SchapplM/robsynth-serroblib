% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRRRRP6
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:11
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRRRP6_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRP6_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_jacobig_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:11:38
	% EndTime: 2019-10-09 23:11:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:11:38
	% EndTime: 2019-10-09 23:11:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:11:38
	% EndTime: 2019-10-09 23:11:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t18 = sin(pkin(6));
	t1 = [0, sin(pkin(12)) * t18, 0, 0, 0, 0; 0, -cos(pkin(12)) * t18, 0, 0, 0, 0; 0, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:11:38
	% EndTime: 2019-10-09 23:11:38
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (6->6), mult. (17->14), div. (0->0), fcn. (28->8), ass. (0->12)
	t64 = sin(pkin(12));
	t66 = sin(pkin(6));
	t74 = t64 * t66;
	t67 = cos(pkin(12));
	t73 = t67 * t66;
	t69 = cos(pkin(6));
	t71 = cos(qJ(2));
	t72 = t69 * t71;
	t70 = sin(qJ(2));
	t68 = cos(pkin(7));
	t65 = sin(pkin(7));
	t1 = [0, t74, -(-t64 * t72 - t67 * t70) * t65 + t68 * t74, 0, 0, 0; 0, -t73, -(-t64 * t70 + t67 * t72) * t65 - t68 * t73, 0, 0, 0; 0, t69, -t66 * t71 * t65 + t69 * t68, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:11:38
	% EndTime: 2019-10-09 23:11:38
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (16->14), mult. (48->33), div. (0->0), fcn. (72->10), ass. (0->18)
	t119 = sin(pkin(12));
	t121 = sin(pkin(6));
	t133 = t119 * t121;
	t120 = sin(pkin(7));
	t132 = t120 * t121;
	t122 = cos(pkin(12));
	t131 = t122 * t121;
	t124 = cos(pkin(6));
	t126 = sin(qJ(2));
	t130 = t124 * t126;
	t128 = cos(qJ(2));
	t129 = t124 * t128;
	t127 = cos(qJ(3));
	t125 = sin(qJ(3));
	t123 = cos(pkin(7));
	t118 = -t119 * t129 - t122 * t126;
	t117 = -t119 * t126 + t122 * t129;
	t1 = [0, t133, -t118 * t120 + t123 * t133, (-t119 * t130 + t122 * t128) * t125 + (-t118 * t123 - t119 * t132) * t127, 0, 0; 0, -t131, -t117 * t120 - t123 * t131, (t119 * t128 + t122 * t130) * t125 + (-t117 * t123 + t120 * t131) * t127, 0, 0; 0, t124, t124 * t123 - t128 * t132, -t124 * t120 * t127 + (-t123 * t127 * t128 + t125 * t126) * t121, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:11:38
	% EndTime: 2019-10-09 23:11:39
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (34->21), mult. (100->47), div. (0->0), fcn. (145->12), ass. (0->29)
	t162 = sin(pkin(12));
	t164 = sin(pkin(6));
	t182 = t162 * t164;
	t163 = sin(pkin(7));
	t181 = t163 * t164;
	t167 = cos(pkin(6));
	t180 = t163 * t167;
	t165 = cos(pkin(12));
	t179 = t165 * t164;
	t166 = cos(pkin(7));
	t173 = cos(qJ(2));
	t178 = t166 * t173;
	t170 = sin(qJ(2));
	t177 = t167 * t170;
	t176 = t167 * t173;
	t158 = -t162 * t170 + t165 * t176;
	t175 = -t158 * t166 + t163 * t179;
	t160 = -t162 * t176 - t165 * t170;
	t174 = t160 * t166 + t162 * t181;
	t172 = cos(qJ(3));
	t171 = cos(qJ(4));
	t169 = sin(qJ(3));
	t168 = sin(qJ(4));
	t161 = -t162 * t177 + t165 * t173;
	t159 = t162 * t173 + t165 * t177;
	t157 = t167 * t166 - t173 * t181;
	t156 = -t160 * t163 + t166 * t182;
	t155 = -t158 * t163 - t166 * t179;
	t1 = [0, t182, t156, t161 * t169 - t174 * t172, (t161 * t172 + t174 * t169) * t168 - t156 * t171, 0; 0, -t179, t155, t159 * t169 + t175 * t172, (t159 * t172 - t175 * t169) * t168 - t155 * t171, 0; 0, t167, t157, -t172 * t180 + (t169 * t170 - t172 * t178) * t164, (t169 * t180 + (t169 * t178 + t170 * t172) * t164) * t168 - t157 * t171, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:11:39
	% EndTime: 2019-10-09 23:11:39
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (34->21), mult. (100->47), div. (0->0), fcn. (145->12), ass. (0->29)
	t197 = sin(pkin(12));
	t199 = sin(pkin(6));
	t217 = t197 * t199;
	t198 = sin(pkin(7));
	t216 = t198 * t199;
	t202 = cos(pkin(6));
	t215 = t198 * t202;
	t200 = cos(pkin(12));
	t214 = t200 * t199;
	t201 = cos(pkin(7));
	t208 = cos(qJ(2));
	t213 = t201 * t208;
	t205 = sin(qJ(2));
	t212 = t202 * t205;
	t211 = t202 * t208;
	t193 = -t197 * t205 + t200 * t211;
	t210 = -t193 * t201 + t198 * t214;
	t195 = -t197 * t211 - t200 * t205;
	t209 = t195 * t201 + t197 * t216;
	t207 = cos(qJ(3));
	t206 = cos(qJ(4));
	t204 = sin(qJ(3));
	t203 = sin(qJ(4));
	t196 = -t197 * t212 + t200 * t208;
	t194 = t197 * t208 + t200 * t212;
	t192 = t202 * t201 - t208 * t216;
	t191 = -t195 * t198 + t201 * t217;
	t190 = -t193 * t198 - t201 * t214;
	t1 = [0, t217, t191, t196 * t204 - t209 * t207, (t196 * t207 + t209 * t204) * t203 - t191 * t206, 0; 0, -t214, t190, t194 * t204 + t210 * t207, (t194 * t207 - t210 * t204) * t203 - t190 * t206, 0; 0, t202, t192, -t207 * t215 + (t204 * t205 - t207 * t213) * t199, (t204 * t215 + (t204 * t213 + t205 * t207) * t199) * t203 - t192 * t206, 0;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end