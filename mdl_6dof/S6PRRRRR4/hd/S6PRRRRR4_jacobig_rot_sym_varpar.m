% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRRRRR4
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
% Datum: 2019-10-09 23:19
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRRRR4_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRR4_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_jacobig_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:19:23
	% EndTime: 2019-10-09 23:19:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:19:23
	% EndTime: 2019-10-09 23:19:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:19:23
	% EndTime: 2019-10-09 23:19:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t18 = sin(pkin(6));
	t1 = [0, sin(pkin(13)) * t18, 0, 0, 0, 0; 0, -cos(pkin(13)) * t18, 0, 0, 0, 0; 0, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:19:23
	% EndTime: 2019-10-09 23:19:23
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
	% StartTime: 2019-10-09 23:19:23
	% EndTime: 2019-10-09 23:19:23
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
	% StartTime: 2019-10-09 23:19:24
	% EndTime: 2019-10-09 23:19:24
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (26->14), mult. (79->33), div. (0->0), fcn. (116->10), ass. (0->21)
	t136 = sin(pkin(13));
	t138 = sin(pkin(6));
	t150 = t136 * t138;
	t137 = sin(pkin(7));
	t149 = t137 * t138;
	t139 = cos(pkin(13));
	t148 = t139 * t138;
	t141 = cos(pkin(6));
	t143 = sin(qJ(2));
	t147 = t141 * t143;
	t145 = cos(qJ(2));
	t146 = t141 * t145;
	t144 = cos(qJ(3));
	t142 = sin(qJ(3));
	t140 = cos(pkin(7));
	t135 = -t136 * t146 - t139 * t143;
	t134 = -t136 * t143 + t139 * t146;
	t133 = -t141 * t137 * t144 + (-t140 * t144 * t145 + t142 * t143) * t138;
	t132 = (-t136 * t147 + t139 * t145) * t142 + (-t135 * t140 - t136 * t149) * t144;
	t131 = (t136 * t145 + t139 * t147) * t142 + (-t134 * t140 + t137 * t148) * t144;
	t1 = [0, t150, -t135 * t137 + t140 * t150, t132, t132, 0; 0, -t148, -t134 * t137 - t140 * t148, t131, t131, 0; 0, t141, t141 * t140 - t145 * t149, t133, t133, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:19:24
	% EndTime: 2019-10-09 23:19:24
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (50->22), mult. (131->47), div. (0->0), fcn. (189->12), ass. (0->33)
	t190 = sin(pkin(13));
	t192 = sin(pkin(6));
	t208 = t190 * t192;
	t191 = sin(pkin(7));
	t207 = t191 * t192;
	t195 = cos(pkin(6));
	t206 = t191 * t195;
	t193 = cos(pkin(13));
	t205 = t193 * t192;
	t194 = cos(pkin(7));
	t199 = cos(qJ(2));
	t204 = t194 * t199;
	t197 = sin(qJ(2));
	t203 = t195 * t197;
	t202 = t195 * t199;
	t183 = -t190 * t197 + t193 * t202;
	t201 = -t183 * t194 + t191 * t205;
	t185 = -t190 * t202 - t193 * t197;
	t200 = t185 * t194 + t190 * t207;
	t198 = cos(qJ(3));
	t196 = sin(qJ(3));
	t189 = qJ(4) + qJ(5);
	t188 = cos(t189);
	t187 = sin(t189);
	t186 = -t190 * t203 + t193 * t199;
	t184 = t190 * t199 + t193 * t203;
	t182 = t195 * t194 - t199 * t207;
	t181 = -t185 * t191 + t194 * t208;
	t180 = -t183 * t191 - t194 * t205;
	t179 = -t198 * t206 + (t196 * t197 - t198 * t204) * t192;
	t178 = t186 * t196 - t200 * t198;
	t177 = t184 * t196 + t201 * t198;
	t1 = [0, t208, t181, t178, t178, (t186 * t198 + t200 * t196) * t187 - t181 * t188; 0, -t205, t180, t177, t177, (t184 * t198 - t201 * t196) * t187 - t180 * t188; 0, t195, t182, t179, t179, (t196 * t206 + (t196 * t204 + t197 * t198) * t192) * t187 - t182 * t188;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end