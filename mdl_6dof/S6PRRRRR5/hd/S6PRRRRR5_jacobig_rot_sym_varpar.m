% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRRRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:21
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRRRR5_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRR5_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_jacobig_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:21:28
	% EndTime: 2019-10-09 23:21:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:21:28
	% EndTime: 2019-10-09 23:21:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:21:28
	% EndTime: 2019-10-09 23:21:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t18 = sin(pkin(6));
	t1 = [0, sin(pkin(13)) * t18, 0, 0, 0, 0; 0, -cos(pkin(13)) * t18, 0, 0, 0, 0; 0, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:21:28
	% EndTime: 2019-10-09 23:21:28
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (6->6), mult. (17->14), div. (0->0), fcn. (28->8), ass. (0->12)
	t64 = sin(pkin(13));
	t66 = sin(pkin(6));
	t74 = t64 * t66;
	t67 = cos(pkin(13));
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
	% StartTime: 2019-10-09 23:21:29
	% EndTime: 2019-10-09 23:21:29
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (16->14), mult. (48->33), div. (0->0), fcn. (72->10), ass. (0->18)
	t119 = sin(pkin(13));
	t121 = sin(pkin(6));
	t133 = t119 * t121;
	t120 = sin(pkin(7));
	t132 = t120 * t121;
	t122 = cos(pkin(13));
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
	% StartTime: 2019-10-09 23:21:29
	% EndTime: 2019-10-09 23:21:29
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (34->21), mult. (100->47), div. (0->0), fcn. (145->12), ass. (0->29)
	t162 = sin(pkin(13));
	t164 = sin(pkin(6));
	t182 = t162 * t164;
	t163 = sin(pkin(7));
	t181 = t163 * t164;
	t167 = cos(pkin(6));
	t180 = t163 * t167;
	t165 = cos(pkin(13));
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
	% StartTime: 2019-10-09 23:21:29
	% EndTime: 2019-10-09 23:21:29
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (52->21), mult. (152->47), div. (0->0), fcn. (218->12), ass. (0->32)
	t180 = sin(pkin(13));
	t182 = sin(pkin(6));
	t200 = t180 * t182;
	t181 = sin(pkin(7));
	t199 = t181 * t182;
	t185 = cos(pkin(6));
	t198 = t181 * t185;
	t183 = cos(pkin(13));
	t197 = t183 * t182;
	t184 = cos(pkin(7));
	t191 = cos(qJ(2));
	t196 = t184 * t191;
	t188 = sin(qJ(2));
	t195 = t185 * t188;
	t194 = t185 * t191;
	t176 = -t180 * t188 + t183 * t194;
	t193 = -t176 * t184 + t181 * t197;
	t178 = -t180 * t194 - t183 * t188;
	t192 = t178 * t184 + t180 * t199;
	t190 = cos(qJ(3));
	t189 = cos(qJ(4));
	t187 = sin(qJ(3));
	t186 = sin(qJ(4));
	t179 = -t180 * t195 + t183 * t191;
	t177 = t180 * t191 + t183 * t195;
	t175 = t185 * t184 - t191 * t199;
	t174 = -t178 * t181 + t184 * t200;
	t173 = -t176 * t181 - t184 * t197;
	t172 = (t187 * t198 + (t187 * t196 + t188 * t190) * t182) * t186 - t175 * t189;
	t171 = (t179 * t190 + t192 * t187) * t186 - t174 * t189;
	t170 = (t177 * t190 - t193 * t187) * t186 - t173 * t189;
	t1 = [0, t200, t174, t179 * t187 - t192 * t190, t171, t171; 0, -t197, t173, t177 * t187 + t193 * t190, t170, t170; 0, t185, t175, -t190 * t198 + (t187 * t188 - t190 * t196) * t182, t172, t172;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end